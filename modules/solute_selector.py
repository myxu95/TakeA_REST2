#!/usr/bin/env python3
"""
REST2 Solute Selector Module - Rewritten based on Pre_analysis.py
Modifies topology files for REST2 enhanced sampling by marking solute atoms
"""

import re
import sys
from pathlib import Path
from typing import List, Set, Dict, Any, Optional
import shutil


class SoluteSelectorError(Exception):
    """Solute selector error"""
    pass


class SoluteSelector:
    """
    Solute selector for REST2 simulations
    Modifies topology files to mark solute atoms for REST2 scaling
    """
    
    def __init__(self, structure_data: Dict[str, Any]):
        """
        Initialize solute selector
        
        Args:
            structure_data: Data from structure analyzer containing:
                - target_atom_indices: Target atom indices
                - nearby_residue_data: List of dicts with chain_id and resid
                - solute_atom_indices: All solute atom indices
                - universe: MDAnalysis Universe object
        """
        self.target_atom_indices = structure_data['target_atom_indices']
        self.nearby_residue_data = structure_data.get('nearby_residue_data', [])
        self.solute_atom_indices = set(structure_data['solute_atom_indices'])
        self.universe = structure_data['universe']
        
        # Get target residue IDs
        self.target_residue_ids = self._get_target_residue_ids()
        
        # All residues involved in REST2 scaling (with chain info)
        self.all_solute_residues = self._get_all_solute_residues()
        
        print(f"✓ SoluteSelector initialization completed:")
        print(f"  Target atoms: {len(self.target_atom_indices)}")
        print(f"  Nearby residues: {len(self.nearby_residue_data)}")
        print(f"  Total solute atoms: {len(self.solute_atom_indices)}")
        print(f"  Solute atom index range: {min(self.solute_atom_indices)} - {max(self.solute_atom_indices)}")
    
    def _get_target_residue_ids(self) -> Set[int]:
        """Get residue IDs for target atoms"""
        target_residues = set()
        
        # If universe is not available, use nearby_residue_data as fallback
        if self.universe is None:
            return {item['resid'] for item in self.nearby_residue_data}
        
        for atom_idx in self.target_atom_indices:
            if atom_idx < len(self.universe.atoms):
                atom = self.universe.atoms[atom_idx]
                if hasattr(atom, 'resid'):
                    target_residues.add(atom.resid)
        return target_residues
    
    def _get_all_solute_residues(self) -> Set[tuple]:
        """Get all solute residues with chain information"""
        all_residues = set()
        
        # Add target residues (assume they're from the target chain)
        for resid in self.target_residue_ids:
            all_residues.add(('target', resid))
        
        # Add nearby residues with chain information
        for item in self.nearby_residue_data:
            chain_id = item['chain_id']
            resid = item['resid']
            all_residues.add((chain_id, resid))
        
        return all_residues
    
    def modify_topology_file(self, input_topology: str, output_topology: str) -> None:
        """
        Modify topology file for REST2 by adding underscores to solute atom types
        Rewritten based on Pre_analysis.py logic
        
        Args:
            input_topology: Input topology file path
            output_topology: Output topology file path
        """
        input_path = Path(input_topology)
        output_path = Path(output_topology)
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input topology file not found: {input_topology}")
        
        print(f"\nStarting topology file modification...")
        print(f"  Input file: {input_topology}")
        print(f"  Output file: {output_topology}")
        
        # Get list of residues to mark
        residues_to_mark = self._get_residues_to_mark()
        print(f"  Residues to mark: {residues_to_mark}")
        
        # Parse and modify topology
        modified_content, modification_stats = self._parse_and_modify_topology(input_path, residues_to_mark)
        
        # Write modified topology
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(modified_content)
        
        # Print modification statistics
        self._print_modification_stats(modification_stats)
        
        # Note: Validation is now handled by the public validate_topology_modification method
        # when called from main.py
    
    def _get_residues_to_mark(self) -> Dict[str, List[int]]:
        """Get mapping of chain IDs to residue IDs that need to be marked"""
        # Create chain-to-residues mapping for accurate topology marking
        chain_residues = {}
        
        # Add target residues to chain C (since target selection is "chainid C")
        target_chain = 'C'  # Target residues are from chain C
        if target_chain not in chain_residues:
            chain_residues[target_chain] = []
        chain_residues[target_chain].extend(sorted(list(self.target_residue_ids)))
        
        # Add nearby residues with their specific chain information
        for item in self.nearby_residue_data:
            chain_id = item['chain_id']
            resid = item['resid']
            
            if chain_id not in chain_residues:
                chain_residues[chain_id] = []
            chain_residues[chain_id].append(resid)
        
        # Sort residue lists for each chain
        for chain_id in chain_residues:
            chain_residues[chain_id].sort()
        
        print(f"    Target residues (Chain C): {sorted(list(self.target_residue_ids))}")
        print(f"    Chain-residue mapping:")
        for chain_id, residues in chain_residues.items():
            print(f"      Chain {chain_id}: {residues}")
        
        return chain_residues
    
    def _parse_and_modify_topology(self, topology_file: Path, residues_to_mark: Dict[str, List[int]]) -> tuple[str, dict]:
        """
        Parse topology file and modify solute atom types
        Rewritten to correctly identify chains and mark only the right residues
        
        Args:
            topology_file: Path to topology file
            residues_to_mark: Mapping of chain IDs to residue IDs to mark
            
        Returns:
            Tuple of (modified topology content, modification statistics)
        """
        # Initialize parsing flags
        found_moleculetype = False
        found_atoms = False
        found_moleculetype_name = False
        found_residue = False
        current_chain_id = None  # Track current chain being parsed
        
        modified_lines = []
        modification_stats = {
            'total_atoms': 0,
            'modified_atoms': 0,
            'molecules_found': set(),
            'molecule_atom_counts': {}
        }
        
        print(f"  Starting topology file parsing...")
        print(f"  Residues to mark by chain: {residues_to_mark}")
        
        with open(topology_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                # Detect [ moleculetype ] section
                if line.strip() == "[ moleculetype ]":
                    found_moleculetype = True
                    found_atoms = False
                    found_moleculetype_name = False
                    found_residue = False
                    current_chain_id = None  # Reset chain ID for new molecule
                    modified_lines.append(line)
                    continue
                
                # Detect [ atoms ] section
                if line.strip() == "[ atoms ]":
                    found_atoms = True
                    modified_lines.append(line)
                    continue
                
                # Detect molecule name and extract chain ID
                if found_moleculetype and not found_atoms:
                    if line.split() and not line.startswith(';'):
                        molecule_name = line.split()[0]
                        found_moleculetype_name = True
                        modification_stats['molecules_found'].add(molecule_name)
                        
                        # Extract chain ID from molecule name
                        current_chain_id = self._extract_chain_id_from_molecule_name(molecule_name)
                        print(f"    Found molecule: {molecule_name} -> Chain: {current_chain_id}")
                        
                    modified_lines.append(line)
                    continue
                
                # Detect residue comment lines
                if found_moleculetype and found_atoms:
                    if line.startswith(";") and "residue" in line:
                        # Parse residue number
                        parts = line.split()
                        if len(parts) >= 3 and parts[2].isdigit():
                            residue_number = int(parts[2])
                            # Reset found_residue flag for this new residue
                            found_residue = False
                            
                            # Check if this residue belongs to the current chain AND needs marking
                            if (current_chain_id and 
                                current_chain_id in residues_to_mark and 
                                residue_number in residues_to_mark[current_chain_id]):
                                found_residue = True
                                print(f"    ✓ Found target residue: {residue_number} (chain {current_chain_id})")
                            else:
                                if current_chain_id:
                                    print(f"    - Skipping residue: {residue_number} (chain {current_chain_id}, not in target list)")
                                else:
                                    print(f"    - Skipping residue: {residue_number} (unknown chain)")
                        
                        modified_lines.append(line)
                        continue
                
                # Detect [ bonds ] section end
                if line.strip() == "[ bonds ]":
                    found_moleculetype = False
                    found_atoms = False
                    found_moleculetype_name = False
                    found_residue = False
                    current_chain_id = None
                    modified_lines.append(line)
                    continue
                
                # Mark atom lines
                if (found_moleculetype and found_atoms and found_moleculetype_name and 
                    not line.startswith(";") and line.strip() and len(line.split()) >= 2):
                    
                    # Only mark atoms if we're currently in a target residue
                    if found_residue:
                        # This is an atom line that needs marking
                        parts = line.split()
                        if parts[0].isdigit():  # Ensure it's an atom line
                            modification_stats['total_atoms'] += 1
                            
                            # Mark atom type (add underscore)
                            atom_type = parts[1]
                            if not atom_type.endswith('_'):
                                # Find position of atom type in the line
                                try:
                                    index = line.index(atom_type) + len(atom_type)
                                    modified_line = line[:index] + "_" + line[index:]
                                    modified_lines.append(modified_line)
                                    modification_stats['modified_atoms'] += 1
                                    print(f"    ✓ Marked atom {parts[0]} (residue {parts[2]} {parts[3]}): {atom_type} → {atom_type}_")
                                except (IndexError, ValueError):
                                    # Fallback method
                                    modified_parts = parts.copy()
                                    modified_parts[1] = atom_type + '_'
                                    modified_line = ' '.join(modified_parts) + '\n'
                                    modified_lines.append(modified_line)
                                    modification_stats['modified_atoms'] += 1
                                    print(f"    ✓ Marked atom {parts[0]} (residue {parts[2]} {parts[3]}): {atom_type} → {atom_type}_ [fallback]")
                            else:
                                # Atom type already has underscore
                                modified_lines.append(modified_line)
                            continue
                    else:
                        # Not in a target residue, just copy the line unchanged
                        modified_lines.append(line)
                        continue
                
                # Other lines remain unchanged
                modified_lines.append(line)
        
        print(f"  Parsing completed:")
        print(f"    Total atoms: {modification_stats['total_atoms']}")
        print(f"    Modified atoms: {modification_stats['modified_atoms']}")
        print(f"    Molecules found: {', '.join(sorted(modification_stats['molecules_found']))}")
        
        return ''.join(modified_lines), modification_stats
    
    def _extract_chain_id_from_molecule_name(self, molecule_name: str) -> str:
        """
        Extract chain ID from GROMACS molecule name
        
        Args:
            molecule_name: Molecule name from topology file
            
        Returns:
            Chain ID (e.g., 'A', 'B', 'C', 'D', 'E') or None if not found
        """
        # Common patterns in GROMACS topology files
        if 'chain' in molecule_name.lower():
            # Extract chain ID from names like "Protein_chain_A", "Chain_B", etc.
            import re
            match = re.search(r'chain[_-]?([A-Z])', molecule_name, re.IGNORECASE)
            if match:
                return match.group(1)
        
        # Handle other naming conventions
        if molecule_name.startswith('Protein_'):
            # Extract from "Protein_chain_A" -> "A"
            parts = molecule_name.split('_')
            if len(parts) >= 3 and parts[1] == 'chain':
                return parts[2]
        
        # If no clear pattern, return None (will skip marking for this molecule)
        return None
    
    def _print_modification_stats(self, stats: dict) -> None:
        """Print topology modification statistics"""
        print(f"\nTopology modification statistics:")
        print(f"  Total atoms: {stats['total_atoms']}")
        print(f"  Modified atoms: {stats['modified_atoms']}")
        if stats['total_atoms'] > 0:
            modification_rate = stats['modified_atoms'] / stats['total_atoms'] * 100
            print(f"  Modification rate: {modification_rate:.1f}%")
        else:
            print(f"  Modification rate: 0%")
    
    def create_rest2_topology_summary(self) -> Dict[str, Any]:
        """
        Create summary of REST2 topology modifications
        
        Returns:
            Dictionary with modification summary
        """
        return {
            'target_residues': list(self.target_residue_ids),
            'nearby_residues': [item['resid'] for item in self.nearby_residue_data],
            'all_solute_residues': list(self.all_solute_residues),
            'target_atoms': len(self.target_atom_indices),
            'solute_atoms': len(self.solute_atom_indices),
            'total_residues': len(self.all_solute_residues)
        }
    
    def print_modification_summary(self) -> None:
        """Print topology modification summary"""
        summary = self.create_rest2_topology_summary()
        
        print("\n" + "="*50)
        print("REST2 Topology Modification Summary")
        print("="*50)
        
        print(f"Target residues: {summary['target_residues']}")
        print(f"Nearby residues: {summary['nearby_residues']}")
        print(f"Total solute residues: {summary['total_residues']}")
        print(f"Target atoms: {summary['target_atoms']}")
        print(f"Solute atoms: {summary['solute_atoms']}")
        
        print("\nResidues to be scaled by REST2:")
        for resid_info in sorted(summary['all_solute_residues']):
            if isinstance(resid_info, tuple):
                chain_id, resid = resid_info
                if resid in summary['target_residues']:
                    print(f"  {chain_id}:{resid} (target)")
                else:
                    print(f"  {chain_id}:{resid} (nearby)")
            else:
                # Handle old format
                resid = resid_info
                if resid in summary['target_residues']:
                    print(f"  {resid} (target)")
                else:
                    print(f"  {resid} (nearby)")
        
        print("="*50)
    
    def validate_topology_modification(self, output_topology: str) -> bool:
        """
        Public method to validate topology modification
        Checks if all selected residues from selected_residues.txt are properly marked
        
        Args:
            output_topology: Path to modified topology file
            
        Returns:
            True if validation successful, False otherwise
        """
        try:
            output_path = Path(output_topology)
            if not output_path.exists():
                print(f"  ⚠ Output topology file not found: {output_topology}")
                return False
            
            print(f"\nValidating topology modification...")
            print(f"  Checking if all selected residues are properly marked...")
            
            # Get the list of residues that should be marked
            residues_to_mark = self._get_residues_to_mark()
            total_expected_residues = sum(len(residues) for residues in residues_to_mark.values())
            
            print(f"  Expected residues to mark: {total_expected_residues}")
            print(f"  Chain breakdown:")
            for chain_id, residues in residues_to_mark.items():
                print(f"    Chain {chain_id}: {len(residues)} residues {residues}")
            
            # Parse the output topology to find marked residues
            marked_residues = set()
            current_residue = None
            
            with open(output_path, 'r') as f:
                for line in f:
                    # Look for residue comment lines
                    if line.startswith(";") and "residue" in line:
                        parts = line.split()
                        if len(parts) >= 3 and parts[2].isdigit():
                            current_residue = int(parts[2])
                    
                    # Look for marked atom lines (with underscore)
                    elif (current_residue is not None and 
                          not line.startswith(";") and 
                          line.strip() and 
                          len(line.split()) >= 2 and
                          line.split()[1].endswith('_')):
                        
                        # This residue has marked atoms
                        marked_residues.add(current_residue)
            
            # Validate marking results
            print(f"\n  Validation results:")
            print(f"    Expected marked residues: {total_expected_residues}")
            print(f"    Actually marked residues: {len(marked_residues)}")
            
            # Check if all expected residues are marked
            all_expected_residues = set()
            for chain_id, residues in residues_to_mark.items():
                all_expected_residues.update(residues)
            
            missing_residues = all_expected_residues - marked_residues
            extra_residues = marked_residues - all_expected_residues
            
            if not missing_residues and not extra_residues:
                print(f"  ✓ All expected residues correctly marked!")
                print(f"  ✓ No extra residues marked!")
                return True
            else:
                if missing_residues:
                    print(f"  ✗ Missing residues: {sorted(missing_residues)}")
                if extra_residues:
                    print(f"  ⚠ Extra residues marked: {sorted(extra_residues)}")
                
                print(f"  ⚠ Marking validation failed")
                return False
                
        except Exception as e:
            print(f"  ⚠ Validation failed: {e}")
            import traceback
            traceback.print_exc()
            return False
    
    def _determine_chain_from_context(self, residue_line: str, resid: int) -> str:
        """
        Determine chain ID from residue context
        This is a heuristic approach since topology files don't always have explicit chain info
        """
        # Try to find chain info in the residue line or nearby context
        # For now, we'll use a simple mapping based on residue ranges
        # This could be improved by parsing molecule names or other context
        
        # Default chain mapping (this is a simplified approach)
        # In practice, you might want to parse this from the topology file structure
        if resid <= 100:  # Assuming first 100 residues are chain A
            return 'A'
        elif resid <= 200:  # Next 100 are chain B
            return 'B'
        elif resid <= 300:  # Next 100 are chain D
            return 'D'
        elif resid <= 400:  # Next 100 are chain E
            return 'E'
        else:
            return 'Unknown'


def main():
    """Test solute selector"""
    try:
        # Create mock structure data for testing
        mock_structure_data = {
            'target_atom_indices': [0, 1, 2, 3],
            'nearby_residue_data': [{'chain_id': 'A', 'resid': 1}, {'chain_id': 'B', 'resid': 2}],
            'solute_atom_indices': [0, 1, 2, 3, 4, 5, 6, 7],
            'universe': None  # Mock universe object
        }
        
        # Test solute selector
        selector = SoluteSelector(mock_structure_data)
        
        print("Testing Solute Selector")
        print("=" * 40)
        
        # Test summary
        selector.print_modification_summary()
        
        print("✓ Solute selector tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()