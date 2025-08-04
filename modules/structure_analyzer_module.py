#!/usr/bin/env python3
"""
REST2 Structure Analysis Module
Analyzes molecular structures and trajectories using MDAnalysis
"""

import MDAnalysis as mda
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Union
import warnings

# Suppress MDAnalysis warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')


class StructureAnalysisError(Exception):
    """Structure analysis error"""
    pass


class StructureAnalyzer:
    """
    Structure analyzer for REST2 simulations
    Handles target identification and nearby residue selection
    """
    
    def __init__(self, structure_file: str, topology_file: str, 
                 trajectory_file: Optional[str] = None):
        """
        Initialize structure analyzer
        
        Args:
            structure_file: Structure file (.gro, .pdb)
            topology_file: Topology file (.top, .tpr)
            trajectory_file: Trajectory file (.xtc, .trr) - optional
        """
        self.structure_file = Path(structure_file)
        self.topology_file = Path(topology_file)
        self.trajectory_file = Path(trajectory_file) if trajectory_file else None
        
        # Validate input files
        self._validate_files()
        
        # Load universe
        self.universe = self._load_universe()
        
        # Analysis results
        self.target_atoms = None
        self.nearby_residues = None
        self.solute_atoms = None
        
    def _validate_files(self) -> None:
        """Validate input files exist"""
        if not self.structure_file.exists():
            raise FileNotFoundError(f"Structure file not found: {self.structure_file}")
        
        if not self.topology_file.exists():
            raise FileNotFoundError(f"Topology file not found: {self.topology_file}")
        
        if self.trajectory_file and not self.trajectory_file.exists():
            raise FileNotFoundError(f"Trajectory file not found: {self.trajectory_file}")
    
    def _load_universe(self) -> mda.Universe:
        """Load MDAnalysis universe"""
        try:
            if self.trajectory_file:
                # Load with trajectory
                universe = mda.Universe(str(self.structure_file), str(self.trajectory_file))
            else:
                # Load structure only
                universe = mda.Universe(str(self.structure_file))
            
            return universe
            
        except Exception as e:
            raise StructureAnalysisError(f"Failed to load structure: {e}")
    
    def identify_target_region(self, target_selection: str) -> mda.AtomGroup:
        """
        Identify target atoms based on selection string
        
        Args:
            target_selection: MDAnalysis selection string (e.g., 'chain A', 'resname LIG')
            
        Returns:
            AtomGroup containing target atoms
        """
        try:
            target_atoms = self.universe.select_atoms(target_selection)
            
            if len(target_atoms) == 0:
                raise StructureAnalysisError(
                    f"No atoms found for selection: {target_selection}"
                )
            
            self.target_atoms = target_atoms
            return target_atoms
            
        except Exception as e:
            raise StructureAnalysisError(f"Target identification failed: {e}")
    
    def find_nearby_residues_static(self, target_atoms: mda.AtomGroup, 
                                  cutoff_distance: float) -> Set[int]:
        """
        Find nearby residues using static structure
        
        Args:
            target_atoms: Target atom group
            cutoff_distance: Cutoff distance in Angstroms
            
        Returns:
            Set of residue IDs within cutoff
        """
        nearby_resids = set()
        
        # Get all protein atoms (excluding target if it's protein)
        protein_atoms = self.universe.select_atoms("protein")
        
        # Calculate distances from target to all protein atoms
        for target_atom in target_atoms:
            target_pos = target_atom.position
            
            for protein_atom in protein_atoms:
                # Skip if same residue (for peptide targets)
                if (hasattr(target_atom, 'resid') and hasattr(protein_atom, 'resid') 
                    and target_atom.resid == protein_atom.resid):
                    continue
                
                distance = np.linalg.norm(target_pos - protein_atom.position)
                
                if distance <= cutoff_distance:
                    nearby_resids.add(protein_atom.resid)
        
        return nearby_resids
    
    def find_nearby_residues_trajectory(self, target_atoms: mda.AtomGroup,
                                      cutoff_distance: float, 
                                      occupancy_threshold: float) -> Set[int]:
        """
        Find nearby residues using trajectory analysis
        
        Args:
            target_atoms: Target atom group
            cutoff_distance: Cutoff distance in Angstroms
            occupancy_threshold: Minimum occupancy fraction (0-1)
            
        Returns:
            Set of residue IDs meeting occupancy threshold
        """
        if not self.trajectory_file:
            raise StructureAnalysisError("Trajectory file required for dynamic analysis")
        
        # Dictionary to track residue contact counts
        residue_contacts = {}
        total_frames = 0
        
        protein_atoms = self.universe.select_atoms("protein")
        
        # Analyze trajectory
        for ts in self.universe.trajectory:
            total_frames += 1
            frame_contacts = set()
            
            # Find contacts in this frame
            for target_atom in target_atoms:
                target_pos = target_atom.position
                
                for protein_atom in protein_atoms:
                    # Skip if same residue
                    if (hasattr(target_atom, 'resid') and hasattr(protein_atom, 'resid')
                        and target_atom.resid == protein_atom.resid):
                        continue
                    
                    distance = np.linalg.norm(target_pos - protein_atom.position)
                    
                    if distance <= cutoff_distance:
                        frame_contacts.add(protein_atom.resid)
            
            # Update contact counts
            for resid in frame_contacts:
                residue_contacts[resid] = residue_contacts.get(resid, 0) + 1
        
        # Filter by occupancy threshold
        nearby_resids = set()
        for resid, count in residue_contacts.items():
            occupancy = count / total_frames
            if occupancy >= occupancy_threshold:
                nearby_resids.add(resid)
        
        return nearby_resids
    
    def analyze_target_and_environment(self, target_selection: str, 
                                     cutoff_distance: float,
                                     use_trajectory: bool = False,
                                     occupancy_threshold: float = 0.5) -> Dict:
        """
        Complete analysis of target region and nearby residues
        
        Args:
            target_selection: Target selection string
            cutoff_distance: Cutoff distance in Angstroms
            use_trajectory: Whether to use trajectory for analysis
            occupancy_threshold: Occupancy threshold for trajectory analysis
            
        Returns:
            Dictionary with analysis results
        """
        # Identify target atoms
        target_atoms = self.identify_target_region(target_selection)
        
        # Find nearby residues
        if use_trajectory:
            nearby_resids = self.find_nearby_residues_trajectory(
                target_atoms, cutoff_distance, occupancy_threshold
            )
        else:
            nearby_resids = self.find_nearby_residues_static(
                target_atoms, cutoff_distance
            )
        
        self.nearby_residues = nearby_resids
        
        # Create solute selection (target + nearby residues)
        solute_atoms = self._create_solute_selection(target_atoms, nearby_resids)
        self.solute_atoms = solute_atoms
        
        # Compile results
        results = {
            'target_atoms': target_atoms,
            'target_atom_count': len(target_atoms),
            'nearby_residues': nearby_resids,
            'nearby_residue_count': len(nearby_resids),
            'solute_atoms': solute_atoms,
            'solute_atom_count': len(solute_atoms),
            'analysis_method': 'trajectory' if use_trajectory else 'static'
        }
        
        if use_trajectory:
            results['total_frames'] = len(self.universe.trajectory)
            results['occupancy_threshold'] = occupancy_threshold
        
        return results
    
    def _create_solute_selection(self, target_atoms: mda.AtomGroup, 
                               nearby_resids: Set[int]) -> mda.AtomGroup:
        """
        Create combined solute selection from target and nearby residues
        
        Args:
            target_atoms: Target atom group
            nearby_resids: Set of nearby residue IDs
            
        Returns:
            Combined atom group for solute
        """
        # Start with target atoms
        solute_selection_parts = []
        
        # Add target selection
        if hasattr(target_atoms, 'residues') and len(target_atoms.residues) > 0:
            # If target is residue-based, include whole residues
            target_resids = [res.resid for res in target_atoms.residues]
            solute_selection_parts.extend([f"resid {resid}" for resid in target_resids])
        else:
            # For non-residue targets (like ligands), we need a different approach
            # This is a simplification - in practice, you might need more sophisticated logic
            solute_selection_parts.append("same residue as (selection_placeholder)")
        
        # Add nearby residues
        if nearby_resids:
            nearby_resids_str = " ".join(map(str, sorted(nearby_resids)))
            solute_selection_parts.append(f"resid {nearby_resids_str}")
        
        # Combine selections
        if len(solute_selection_parts) > 1:
            combined_selection = " or ".join([f"({part})" for part in solute_selection_parts])
        else:
            combined_selection = solute_selection_parts[0]
        
        # Handle the placeholder case for non-residue targets
        if "selection_placeholder" in combined_selection:
            # Get atom indices for target atoms and create selection
            target_indices = " ".join(map(str, target_atoms.indices))
            combined_selection = combined_selection.replace(
                "same residue as (selection_placeholder)", 
                f"index {target_indices}"
            )
        
        try:
            solute_atoms = self.universe.select_atoms(combined_selection)
            return solute_atoms
        except Exception as e:
            # Fallback: just combine target atoms with nearby residue atoms
            nearby_atoms = self.universe.select_atoms(f"resid {' '.join(map(str, nearby_resids))}")
            return target_atoms + nearby_atoms
    
    def get_selected_residues_info(self, analysis_results: Optional[Dict] = None) -> List[Dict]:
        """
        Get detailed information about selected residues for user review
        
        Args:
            analysis_results: Analysis results (uses stored results if None)
            
        Returns:
            List of dictionaries with residue information
        """
        if analysis_results is None:
            if self.nearby_residues is None:
                raise StructureAnalysisError("No analysis results available")
            nearby_residues = self.nearby_residues
        else:
            nearby_residues = analysis_results['nearby_residues']
        
        residue_info = []
        
        for resid in sorted(nearby_residues):
            try:
                residue = self.universe.select_atoms(f"resid {resid}").residues[0]
                info = {
                    'resid': resid,
                    'resname': residue.resname,
                    'chain': getattr(residue, 'chainID', 'Unknown'),
                    'atom_count': len(residue.atoms)
                }
                residue_info.append(info)
            except (IndexError, AttributeError):
                # Handle cases where residue might not be found
                residue_info.append({
                    'resid': resid,
                    'resname': 'Unknown',
                    'chain': 'Unknown',
                    'atom_count': 0
                })
        
        return residue_info
    
    def get_solute_selection_data(self, analysis_results: Optional[Dict] = None) -> Dict:
        """
        Get solute selection data for next module
        
        Args:
            analysis_results: Analysis results (uses stored results if None)
            
        Returns:
            Dictionary with selection data for solute selector module
        """
        if analysis_results is None:
            if self.solute_atoms is None or self.target_atoms is None:
                raise StructureAnalysisError("No analysis results available")
            target_atoms = self.target_atoms
            nearby_residues = self.nearby_residues
            solute_atoms = self.solute_atoms
        else:
            target_atoms = analysis_results['target_atoms']
            nearby_residues = analysis_results['nearby_residues']
            solute_atoms = analysis_results['solute_atoms']
        
        return {
            'target_atom_indices': target_atoms.indices,
            'nearby_residue_ids': list(nearby_residues),
            'solute_atom_indices': solute_atoms.indices,
            'universe': self.universe
        }
    
    def print_selected_residues(self, analysis_results: Optional[Dict] = None) -> None:
        """
        Print selected residues for user review
        
        Args:
            analysis_results: Analysis results (uses stored results if None)
        """
        residue_info = self.get_selected_residues_info(analysis_results)
        
        print("\n" + "="*60)
        print("Selected Residues for REST2 Enhanced Sampling")
        print("="*60)
        print(f"{'ResID':<8} {'ResName':<8} {'Chain':<8} {'Atoms':<8}")
        print("-"*60)
        
        for info in residue_info:
            print(f"{info['resid']:<8} {info['resname']:<8} "
                  f"{info['chain']:<8} {info['atom_count']:<8}")
        
        print("-"*60)
        print(f"Total selected residues: {len(residue_info)}")
        
        if analysis_results:
            print(f"Target atoms: {analysis_results['target_atom_count']}")
            print(f"Total solute atoms: {analysis_results['solute_atom_count']}")
            print(f"Analysis method: {analysis_results['analysis_method']}")
        
        print("="*60 + "\n")


def main():
    """Test structure analyzer"""
    try:
        # Example usage
        structure_file = "example/MD_results/md.gro"
        topology_file = "example/MD_results/topol.top"
        trajectory_file = "example/MD_results/md.xtc"  # Optional
        
        # Initialize analyzer
        analyzer = StructureAnalyzer(structure_file, topology_file, trajectory_file)
        
        # Analyze peptide target
        results = analyzer.analyze_target_and_environment(
            target_selection="chain A",
            cutoff_distance=6.0,
            use_trajectory=False
        )
        
        # Print selected residues for user review
        analyzer.print_selected_residues(results)
        
        # Get data for next module
        solute_data = analyzer.get_solute_selection_data(results)
        print(f"Data prepared for solute selector:")
        print(f"  - Target atoms: {len(solute_data['target_atom_indices'])}")
        print(f"  - Nearby residues: {len(solute_data['nearby_residue_ids'])}")
        print(f"  - Total solute atoms: {len(solute_data['solute_atom_indices'])}")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()