#!/usr/bin/env python3
"""
REST2 Diagnostic Tool
Check for potential issues in REST2 setup that could cause large energy differences
"""

import re
import sys
from pathlib import Path


def check_marked_atoms(topology_file):
    """Check the number and distribution of marked atoms"""
    print("="*60)
    print("Checking Marked Atoms")
    print("="*60)

    total_atoms = 0
    marked_atoms = 0
    marked_by_molecule = {}
    current_molecule = None
    in_atoms = False

    with open(topology_file, 'r') as f:
        for line in f:
            stripped = line.strip()

            # Track molecule
            if stripped.startswith('[ moleculetype ]'):
                in_moleculetype = True
                current_molecule = None
                continue

            if 'in_moleculetype' in locals() and in_moleculetype:
                if not stripped.startswith(';') and not stripped.startswith('[') and stripped:
                    parts = stripped.split()
                    if parts:
                        current_molecule = parts[0]
                        marked_by_molecule[current_molecule] = {'total': 0, 'marked': 0}
                        in_moleculetype = False

            # Track atoms
            if stripped.startswith('[ atoms ]'):
                in_atoms = True
                continue

            if in_atoms and stripped.startswith('['):
                in_atoms = False

            if in_atoms and not stripped.startswith(';') and stripped:
                parts = stripped.split()
                if len(parts) >= 8:
                    total_atoms += 1
                    atom_type = parts[1]

                    if current_molecule:
                        marked_by_molecule[current_molecule]['total'] += 1

                    if atom_type.endswith('_'):
                        marked_atoms += 1
                        if current_molecule:
                            marked_by_molecule[current_molecule]['marked'] += 1

    print(f"Total atoms: {total_atoms}")
    print(f"Marked atoms: {marked_atoms}")
    print(f"Marking rate: {marked_atoms/total_atoms*100:.2f}%")
    print(f"\nPer-molecule breakdown:")
    for mol, stats in marked_by_molecule.items():
        if stats['marked'] > 0:
            rate = stats['marked']/stats['total']*100 if stats['total'] > 0 else 0
            print(f"  {mol}: {stats['marked']}/{stats['total']} ({rate:.1f}%)")
    print()


def check_charge_scaling(topology_marked, topology_scaled, lambda_val):
    """Check if charges are properly scaled"""
    print("="*60)
    print(f"Checking Charge Scaling (λ={lambda_val})")
    print("="*60)

    # Extract marked atom charges from marked topology
    marked_charges = {}
    with open(topology_marked, 'r') as f:
        in_atoms = False
        atom_idx = 0
        for line in f:
            stripped = line.strip()
            if stripped.startswith('[ atoms ]'):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith('['):
                in_atoms = False
            if in_atoms and not stripped.startswith(';') and stripped:
                parts = stripped.split()
                if len(parts) >= 8:
                    atom_type = parts[1]
                    if atom_type.endswith('_'):
                        charge = float(parts[6])
                        marked_charges[atom_idx] = charge
                    atom_idx += 1

    # Extract scaled charges
    scaled_charges = {}
    with open(topology_scaled, 'r') as f:
        in_atoms = False
        atom_idx = 0
        for line in f:
            stripped = line.strip()
            if stripped.startswith('[ atoms ]'):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith('['):
                in_atoms = False
            if in_atoms and not stripped.startswith(';') and stripped:
                parts = stripped.split()
                if len(parts) >= 7:
                    charge = float(parts[6])
                    if atom_idx in marked_charges:
                        scaled_charges[atom_idx] = charge
                    atom_idx += 1

    # Check scaling
    import math
    expected_factor = math.sqrt(lambda_val)
    errors = []

    for idx in list(marked_charges.keys())[:10]:  # Check first 10 marked atoms
        if idx in scaled_charges:
            orig = marked_charges[idx]
            scaled = scaled_charges[idx]
            expected = orig * expected_factor
            diff = abs(scaled - expected)

            if diff > 1e-5:
                errors.append(f"  Atom {idx}: {orig:.6f} → {scaled:.6f} (expected {expected:.6f}, diff {diff:.6e})")

    if errors:
        print("⚠ Charge scaling errors found:")
        for err in errors:
            print(err)
    else:
        print(f"✓ Charge scaling correct (factor = √{lambda_val} = {expected_factor:.6f})")
        print(f"  Checked {min(10, len(marked_charges))} marked atoms")
    print()


def check_temperature_ladder(output_dir, n_replicas):
    """Check temperature ladder for reasonable overlaps"""
    print("="*60)
    print("Checking Temperature Ladder")
    print("="*60)

    # Read temperatures from replica_info files
    temps = []
    for i in range(n_replicas):
        info_file = Path(output_dir) / f"replica_{i}" / "replica_info.txt"
        if info_file.exists():
            with open(info_file, 'r') as f:
                for line in f:
                    if 'Temperature:' in line:
                        temp = float(line.split(':')[1].strip().split()[0])
                        temps.append((i, temp))
                        break

    if len(temps) < 2:
        print("⚠ Could not read temperature ladder")
        return

    temps.sort(key=lambda x: x[1])

    print(f"Temperature range: {temps[0][1]:.1f}K - {temps[-1][1]:.1f}K")
    print(f"Number of replicas: {len(temps)}")

    # Check spacing
    spacings = []
    for i in range(len(temps)-1):
        spacing = temps[i+1][1] - temps[i][1]
        spacings.append(spacing)

    avg_spacing = sum(spacings) / len(spacings)
    max_spacing = max(spacings)
    min_spacing = min(spacings)

    print(f"Average spacing: {avg_spacing:.2f}K")
    print(f"Max spacing: {max_spacing:.2f}K")
    print(f"Min spacing: {min_spacing:.2f}K")

    # Warning if spacing is too large
    if max_spacing > 15:
        print(f"⚠ WARNING: Large temperature spacing ({max_spacing:.1f}K) may cause poor exchange rates")
    else:
        print("✓ Temperature spacing looks reasonable")
    print()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python diagnose_rest2.py <output_dir>")
        print("Example: python diagnose_rest2.py rest2_samples/1bd2/02_rest2_300-450K_24rep")
        sys.exit(1)

    output_dir = Path(sys.argv[1])

    # Check marked topology
    marked_top = output_dir / "rest2_topol_marked.top"
    if marked_top.exists():
        check_marked_atoms(marked_top)
    else:
        print(f"⚠ Marked topology not found: {marked_top}")

    # Check charge scaling for replica 10 (mid-range)
    replica_10_top = output_dir / "replica_10" / "input" / "topol.top"
    replica_10_info = output_dir / "replica_10" / "replica_info.txt"

    if marked_top.exists() and replica_10_top.exists() and replica_10_info.exists():
        # Read lambda from replica_info
        lambda_val = None
        with open(replica_10_info, 'r') as f:
            for line in f:
                if 'Scaling factor:' in line:
                    lambda_val = float(line.split(':')[1].strip())
                    break

        if lambda_val:
            check_charge_scaling(marked_top, replica_10_top, lambda_val)

    # Check temperature ladder
    check_temperature_ladder(output_dir, 24)

    print("="*60)
    print("Diagnostic Complete")
    print("="*60)
