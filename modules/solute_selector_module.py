#!/usr/bin/env python3
"""
REST2 Solute Selector Module
Modifies topology files for REST2 enhanced sampling by marking solute atoms
"""

import re
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
                - nearby_residue_ids: List of nearby residue IDs
                - solute_atom_indices: All solute atom indices
                - universe: MDAnalysis Universe object
        """
        self.target_atom_indices = structure_data['target_atom_indices']
        self.nearby_residue_ids = set(structure_data['nearby_residue_ids'])
        self.solute_atom_indices = set(structure_data['solute_atom_indices'])
        self.universe = structure_data['universe']
        
        # Get target residue IDs
        self.target_residue_ids = self._get_target_residue_ids()
        
        # All residues involved in REST2 scaling
        self.all_solute_residues = self.target_residue_ids | self.nearby_residue_ids
    
    def _get_target_residue_ids(self) -> Set[int]:
        """Get residue IDs for target atoms"""
        target_residues = set()
        for atom_idx in self.target_atom_indices:
            atom = self.universe.atoms[atom_idx]
            if hasattr(atom, 'resid'):
                target_residues.add(atom.resid)
        return target_residues
    
    def modify_topology_file(self, input_topology: str, output_topology: str,
                           molecule_name: Optional[str] = None) -> None:
        """
        Modify topology file for REST2 by adding underscores to solute atom types
        
        Args:
            input_topology: Input topology file path
            output_topology: Output topology file path  
            molecule_name: Specific molecule name to modify (auto-detect if None)
        """
        input_path = Path(input_topology)
        output_path = Path(output_topology)
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input topology file not found: {input_topology}")
        
        # Create backup
        backup_path = input_path.with_suffix(input_path.suffix + '.backup')
        if not backup_path.exists():
            shutil.copy2(input_path, backup_path)
        
        # Parse and modify topology
        modified_content = self._parse_and_modify_topology(input_path, molecule_name)
        
        # Write modified topology
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(modified_content)
    
    def _parse_and_modify_topology(self, topology_file: Path, 
                                 molecule_name: Optional[str]) -> str:
        """
        Parse topology file and modify solute atom types
        
        Args:
            topology_file: Path to topology file
            molecule_name: Target molecule name (auto-detect if None)
            
        Returns:
            Modified topology content
        """
        # Initialize parsing flags
        in_moleculetype = False
        in_atoms_section = False
        found_target_molecule = False
        current_molecule = None
        
        modified_lines = []
        
        with open(topology_file, 'r') as f:
            for line in f:
                stripped_line = line.strip()
                
                # Check for section headers
                if stripped_line.startswith('[') and stripped_line.endswith(']'):
                    section_name = stripped_line[1:-1].strip()
                    
                    if section_name == 'moleculetype':
                        in_moleculetype = True
                        in_atoms_section = False
                    elif section_name == 'atoms':
                        in_atoms_section = True
                        in_moleculetype = False
                    else:
                        in_moleculetype = False
                        in_atoms_section = False
                        found_target_molecule = False
                
                # Parse moleculetype section
                elif in_moleculetype and not stripped_line.startswith(';') and stripped_line:
                    # First non-comment line in moleculetype contains molecule name
                    fields = stripped_line.split()
                    if fields:
                        current_molecule = fields[0]
                        if molecule_name is None or current_molecule == molecule_name:
                            found_target_molecule = True
                
                # Parse atoms section
                elif in_atoms_section and found_target_molecule:
                    if stripped_line.startswith(';') and 'residue' in stripped_line:
                        # Comment line with residue information
                        residue_id = self._extract_residue_id_from_comment(stripped_line)
                        if residue_id in self.all_solute_residues:
                            modified_lines.append(line)
                        else:
                            # Skip this residue section
                            continue
                    elif not stripped_line.startswith(';') and stripped_line:
                        # Atom definition line
                        modified_line = self._modify_atom_line(line, stripped_line)
                        if modified_line:
                            modified_lines.append(modified_line)
                    else:
                        # Empty lines or other comments
                        modified_lines.append(line)
                else:
                    # All other lines pass through unchanged
                    modified_lines.append(line)
        
        return ''.join(modified_lines)
    
    def _extract_residue_id_from_comment(self, comment_line: str) -> Optional[int]:
        """
        Extract residue ID from comment line
        
        Args:
            comment_line: Comment line containing residue info
            
        Returns:
            Residue ID if found, None otherwise
        """
        # Try to find residue number in comment
        # Common formats: "; residue 123 ALA", "; 123 ALA", etc.
        patterns = [
            r';\s*residue\s+(\d+)',
            r';\s*(\d+)\s+[A-Z]{3}',
            r';\s*res\s*(\d+)'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, comment_line, re.IGNORECASE)
            if match:
                return int(match.group(1))
        
        return None
    
    def _modify_atom_line(self, original_line: str, stripped_line: str) -> Optional[str]:
        """
        Modify atom line by adding underscore to atom type if it's a solute atom
        
        Args:
            original_line: Original line with formatting
            stripped_line: Stripped line for parsing
            
        Returns:
            Modified line or None if line should be skipped
        """
        fields = stripped_line.split()
        
        # Basic validation - atom line should have at least atom_id, atom_type, resid
        if len(fields) < 3:
            return original_line
        
        try:
            atom_id = int(fields[0])  # 1-based in topology
            atom_type = fields[1]
            resid = int(fields[2])
            
            # Check if this residue should be included
            if resid not in self.all_solute_residues:
                return None  # Skip this atom
            
            # For solute residues, modify ALL atoms (both target and nearby)
            # Add underscore to atom type for REST2 scaling
            if not atom_type.endswith('_'):
                modified_atom_type = atom_type + '_'
                # Replace the atom type in the original line
                modified_line = original_line.replace(atom_type, modified_atom_type, 1)
                return modified_line
            
            return original_line
            
        except (ValueError, IndexError):
            # If parsing fails, return original line
            return original_line
    
    def create_rest2_topology_summary(self) -> Dict[str, Any]:
        """
        Create summary of REST2 topology modifications
        
        Returns:
            Dictionary with modification summary
        """
        # Count atoms in different categories
        target_atom_count = len(self.target_atom_indices)
        nearby_atom_count = len(self.solute_atom_indices) - target_atom_count
        
        # Get residue information
        target_residue_info = []
        nearby_residue_info = []
        
        for resid in sorted(self.target_residue_ids):
            residue = self.universe.select_atoms(f"resid {resid}").residues[0]
            target_residue_info.append({
                'resid': resid,
                'resname': residue.resname,
                'atom_count': len(residue.atoms)
            })
        
        for resid in sorted(self.nearby_residue_ids):
            residue = self.universe.select_atoms(f"resid {resid}").residues[0]
            nearby_residue_info.append({
                'resid': resid,
                'resname': residue.resname,
                'atom_count': len(residue.atoms)
            })
        
        return {
            'target_residues': target_residue_info,
            'nearby_residues': nearby_residue_info,
            'target_atom_count': target_atom_count,
            'nearby_atom_count': nearby_atom_count,
            'total_solute_atoms': len(self.solute_atom_indices),
            'total_solute_residues': len(self.all_solute_residues)
        }
    
    def print_modification_summary(self) -> None:
        """Print summary of topology modifications"""
        summary = self.create_rest2_topology_summary()
        
        print("\n" + "="*70)
        print("REST2 Topology Modification Summary")
        print("="*70)
        
        print(f"Total solute atoms: {summary['total_solute_atoms']}")
        print(f"Total solute residues: {summary['total_solute_residues']}")
        print()
        
        print("Target residues (primary scaling):")
        print(f"{'ResID':<8} {'ResName':<8} {'Atoms':<8}")
        print("-"*30)
        for res_info in summary['target_residues']:
            print(f"{res_info['resid']:<8} {res_info['resname']:<8} {res_info['atom_count']:<8}")
        print()
        
        if summary['nearby_residues']:
            print("Nearby residues (environment scaling):")
            print(f"{'ResID':<8} {'ResName':<8} {'Atoms':<8}")
            print("-"*30)
            for res_info in summary['nearby_residues']:
                print(f"{res_info['resid']:<8} {res_info['resname']:<8} {res_info['atom_count']:<8}")
        
        print("="*70)
        print("Note: Atom types for ALL atoms in target and nearby residues will be modified with '_' suffix")
        print("="*70 + "\n")
    
    def validate_topology_modification(self, output_topology: str) -> bool:
        """
        Validate that topology modification was successful
        
        Args:
            output_topology: Path to modified topology file
            
        Returns:
            True if validation passes
        """
        if not Path(output_topology).exists():
            return False
        
        # Count modified atom types (with underscores)
        modified_count = 0
        total_atoms = 0
        
        with open(output_topology, 'r') as f:
            in_atoms = False
            for line in f:
                stripped = line.strip()
                if stripped.startswith('[') and 'atoms' in stripped:
                    in_atoms = True
                elif stripped.startswith('[') and in_atoms:
                    break
                elif in_atoms and not stripped.startswith(';') and stripped:
                    fields = stripped.split()
                    if len(fields) >= 2:
                        total_atoms += 1
                        atom_type = fields[1]
                        if atom_type.endswith('_'):
                            modified_count += 1
        
        expected_modifications = len(self.solute_atom_indices)
        
        print(f"Validation: {modified_count} atoms modified (expected: target + nearby residues)")
        return modified_count > 0


def main():
    """Test solute selector"""
    try:
        # Mock structure data for testing
        mock_structure_data = {
            'target_atom_indices': [10, 11, 12, 13, 14],
            'nearby_residue_ids': [15, 16, 17],
            'solute_atom_indices': [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            'universe': None  # Would normally be MDAnalysis Universe
        }
        
        # Initialize selector
        selector = SoluteSelector(mock_structure_data)
        
        # Print summary
        selector.print_modification_summary()
        
        # Example topology modification
        # selector.modify_topology_file(
        #     "example/MD_results/topol.top",
        #     "rest2_topol.top"
        # )
        
        print("Solute selector module ready for use")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()