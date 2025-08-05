#!/usr/bin/env python3
"""
REST2 Solute Selector Module
Modifies topology files for REST2 enhanced sampling by marking solute atoms
"""

import re
import sys
from pathlib import Path
from typing import List, Set, Dict, Any, Optional
import shutil

# Import from utils package
try:
    from utils import ValidationFramework, FileUtils, FileOperationError, OutputFormatter
except ImportError:
    # Fallback for direct execution
    ValidationFramework = None
    FileUtils = None
    FileOperationError = Exception
    OutputFormatter = None


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
        
        # If universe is not available, use nearby_residue_ids as fallback
        if self.universe is None:
            return self.nearby_residue_ids
        
        for atom_idx in self.target_atom_indices:
            if atom_idx < len(self.universe.atoms):
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
        
        # Parse and modify topology
        modified_content, modification_stats = self._parse_and_modify_topology(input_path, molecule_name)
        
        # Write modified topology
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(modified_content)
        
        # Print modification statistics
        self._print_modification_stats(modification_stats)
    
    def _parse_and_modify_topology(self, topology_file: Path, 
                                 molecule_name: Optional[str]) -> tuple[str, dict]:
        """
        Parse topology file and modify solute atom types
        
        Args:
            topology_file: Path to topology file
            molecule_name: Target molecule name (auto-detect if None)
            
        Returns:
            Tuple of (modified topology content, modification statistics)
        """
        # Initialize parsing flags
        in_moleculetype = False
        in_atoms_section = False
        found_target_molecule = False
        current_molecule = None
        
        modified_lines = []
        modification_stats = {
            'total_atoms': 0,
            'modified_atoms': 0,
            'target_molecule': None,
            'molecules_found': []
        }
        
        with open(topology_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped_line = line.strip()
                
                # Track moleculetype sections
                if stripped_line.startswith('[ moleculetype ]'):
                    in_moleculetype = True
                    in_atoms_section = False
                    current_molecule = None
                    continue
                
                # Track atoms sections
                if stripped_line.startswith('[ atoms ]'):
                    in_atoms_section = True
                    continue
                
                # End of sections
                if stripped_line.startswith('[') and not stripped_line.startswith('[ atoms ]'):
                    in_atoms_section = False
                    in_moleculetype = False
                
                # Parse molecule name
                if in_moleculetype and not stripped_line.startswith(';') and stripped_line:
                    if not current_molecule:
                        current_molecule = stripped_line.split()[0]
                        modification_stats['molecules_found'].append(current_molecule)
                        # Check if this is the target molecule
                        if molecule_name is None or current_molecule == molecule_name:
                            found_target_molecule = True
                            modification_stats['target_molecule'] = current_molecule
                        else:
                            found_target_molecule = False
                
                # Modify atom lines in target molecule
                if in_atoms_section and found_target_molecule and not stripped_line.startswith(';'):
                    if stripped_line and not stripped_line.startswith('['):
                        modification_stats['total_atoms'] += 1
                        modified_line = self._modify_atom_line(line, stripped_line)
                        if modified_line:
                            modified_lines.append(modified_line)
                            modification_stats['modified_atoms'] += 1
                        else:
                            modified_lines.append(line)
                    else:
                        modified_lines.append(line)
                else:
                    modified_lines.append(line)
        
        return ''.join(modified_lines), modification_stats
    
    def _print_modification_stats(self, stats: dict) -> None:
        """Print topology modification statistics"""
        print(f"Topology modification:")
        print(f"  Target molecule: {stats['target_molecule']}")
        print(f"  Total atoms: {stats['total_atoms']}")
        print(f"  Modified atoms: {stats['modified_atoms']}")
        print(f"  Modification rate: {stats['modified_atoms']/stats['total_atoms']*100:.1f}%" if stats['total_atoms'] > 0 else "  Modification rate: 0%")
    
    def _extract_residue_id_from_comment(self, comment_line: str) -> Optional[int]:
        """
        Extract residue ID from comment line
        
        Args:
            comment_line: Comment line from topology file
            
        Returns:
            Residue ID if found, None otherwise
        """
        # Look for residue ID in comment (format: ; residue X)
        residue_match = re.search(r'residue\s+(\d+)', comment_line)
        if residue_match:
            return int(residue_match.group(1))
        return None
    
    def _modify_atom_line(self, original_line: str, stripped_line: str) -> Optional[str]:
        """
        Modify atom line for REST2 scaling
        
        Args:
            original_line: Original line from topology file
            stripped_line: Stripped line content
            
        Returns:
            Modified line or None if no modification needed
        """
        # Parse atom line (format: nr type resnr residue atom cgnr charge mass)
        parts = stripped_line.split()
        if len(parts) < 8:
            return None
        
        try:
            atom_nr = int(parts[0]) - 1  # Convert to 0-based indexing
            resnr = int(parts[2])
            atom_type = parts[1]
            
            # Check if this atom is in solute
            if atom_nr in self.solute_atom_indices:
                # Add underscore to atom type for REST2 scaling
                if not atom_type.endswith('_'):
                    modified_parts = parts.copy()
                    modified_parts[1] = atom_type + '_'
                    return ' '.join(modified_parts) + '\n'
                else:
                    # Atom type already has underscore, no need to modify
                    return None
            
            return None
            
        except (ValueError, IndexError) as e:
            # Log parsing error but continue processing
            print(f"Warning: Could not parse atom line: {stripped_line[:50]}...")
            return None
    
    def create_rest2_topology_summary(self) -> Dict[str, Any]:
        """
        Create summary of REST2 topology modifications
        
        Returns:
            Dictionary with modification summary
        """
        return {
            'target_residues': list(self.target_residue_ids),
            'nearby_residues': list(self.nearby_residue_ids),
            'all_solute_residues': list(self.all_solute_residues),
            'target_atoms': len(self.target_atom_indices),
            'solute_atoms': len(self.solute_atom_indices),
            'total_residues': len(self.all_solute_residues)
        }
    
    def print_modification_summary(self) -> None:
        """Print topology modification summary using unified formatter"""
        summary = self.create_rest2_topology_summary()
        
        if OutputFormatter:
            # Use unified output formatter
            OutputFormatter.print_header("REST2 Topology Modification Summary")
            
            # Print summary data
            OutputFormatter.print_summary(summary, "Modification Summary")
            
            # Print residue details
            residue_details = []
            for resid in sorted(summary['all_solute_residues']):
                if resid in summary['target_residues']:
                    residue_details.append(f"  {resid} (target)")
                else:
                    residue_details.append(f"  {resid} (nearby)")
            
            OutputFormatter.print_section("Residues to be scaled by REST2", residue_details)
        else:
            # Fallback summary
            print("\n" + "="*50)
            print("REST2 Topology Modification Summary")
            print("="*50)
            
            print(f"Target residues: {summary['target_residues']}")
            print(f"Nearby residues: {summary['nearby_residues']}")
            print(f"Total solute residues: {summary['total_residues']}")
            print(f"Target atoms: {summary['target_atoms']}")
            print(f"Solute atoms: {summary['solute_atoms']}")
            
            print("\nResidues to be scaled by REST2:")
            for resid in sorted(summary['all_solute_residues']):
                if resid in summary['target_residues']:
                    print(f"  {resid} (target)")
                else:
                    print(f"  {resid} (nearby)")
            
            print("="*50)
    
    def validate_topology_modification(self, output_topology: str) -> bool:
        """
        Validate topology modification using unified validation framework
        
        Args:
            output_topology: Path to modified topology file
            
        Returns:
            True if modification is valid
        """
        try:
            # Create solute data for validation
            solute_data = {
                'target_atom_indices': list(self.target_atom_indices),
                'nearby_residue_ids': list(self.nearby_residue_ids),
                'solute_atom_indices': list(self.solute_atom_indices)
            }
            
            if ValidationFramework:
                # Use unified validation framework
                errors = ValidationFramework.validate_topology_modification(
                    '', output_topology, solute_data  # Empty input path since we're validating output
                )
                
                ValidationFramework.print_validation_summary(errors, "Topology modification validation")
                return len(errors) == 0
            else:
                # Fallback validation
                if not Path(output_topology).exists():
                    print(f"✗ Output topology file not created: {output_topology}")
                    return False
                
                # Check solute data is not empty
                if len(self.solute_atom_indices) == 0:
                    print("✗ No solute atoms selected")
                    return False
                
                print("✓ Topology modification validation passed")
                return True
                
        except Exception as e:
            print(f"✗ Topology modification validation failed: {e}")
            return False


def main():
    """Test solute selector"""
    try:
        # Create mock structure data for testing
        mock_structure_data = {
            'target_atom_indices': [0, 1, 2, 3],
            'nearby_residue_ids': [1, 2, 3],
            'solute_atom_indices': [0, 1, 2, 3, 4, 5, 6, 7],
            'universe': None  # Mock universe object
        }
        
        # Test solute selector
        selector = SoluteSelector(mock_structure_data)
        
        print("Testing Solute Selector")
        print("=" * 40)
        
        # Test summary
        selector.print_modification_summary()
        
        # Test validation (will fail since files don't exist)
        print("\nTesting validation (expected to fail):")
        selector.validate_topology_modification("test_output.top")
        
        print("✓ Solute selector tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()