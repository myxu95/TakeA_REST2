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
import datetime

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
            # Priority 1: If we have TPR file, try to create PDB with chain information
            if self.topology_file.suffix == '.tpr':
                pdb_file = self.structure_file.with_suffix('.pdb')
                if not pdb_file.exists():
                    self._extract_pdb_with_chains()
                
                # Load PDB file if it exists
                if pdb_file.exists():
                    if self.trajectory_file:
                        universe = mda.Universe(str(pdb_file), str(self.trajectory_file))
                    else:
                        universe = mda.Universe(str(pdb_file))
                    return universe
            
            # Fallback to original method
            if self.trajectory_file:
                # Load with trajectory: GRO + XTC
                universe = mda.Universe(str(self.structure_file), str(self.trajectory_file))
            else:
                # Load structure only: GRO file
                universe = mda.Universe(str(self.structure_file))
            
            return universe
            
        except Exception as e:
            raise StructureAnalysisError(f"Failed to load structure: {e}")
    
    def _extract_pdb_with_chains(self) -> None:
        """Extract PDB file with chain information from GRO and TPR"""
        try:
            import subprocess
            
            pdb_file = self.structure_file.with_suffix('.pdb')
            
            # Run gmx trjconv to convert GRO to PDB with chain information
            cmd = [
                'gmx', 'trjconv',
                '-f', str(self.structure_file.absolute()),  # Input GRO file with absolute path
                '-s', str(self.topology_file.absolute()),   # Input TPR file with absolute path
                '-o', str(pdb_file.absolute()),            # Output PDB file with absolute path
                '-pbc', 'mol'
            ]
            
            print(f"Converting GRO to PDB with chain information...")
            print(f"Input GRO: {self.structure_file.absolute()}")
            print(f"Input TPR: {self.topology_file.absolute()}")
            print(f"Output PDB: {pdb_file.absolute()}")
            
            # Run command with input "0" (System)
            result = subprocess.run(
                cmd,
                input="0\n",
                capture_output=True,
                text=True,
                cwd=self.structure_file.parent
            )
            
            if result.returncode == 0 and pdb_file.exists():
                print(f"✓ Successfully converted GRO to PDB: {pdb_file}")
            else:
                print(f"Warning: Could not convert GRO to PDB")
                if result.stderr:
                    print(f"Error: {result.stderr}")
                if result.stdout:
                    print(f"Output: {result.stdout}")
            
        except Exception as e:
            print(f"Warning: Could not convert GRO to PDB: {e}")
            return
    
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
                                   cutoff_distance: float) -> Set[tuple]:
        """
        Find nearby residues using improved atom-to-atom distance calculation
        
        Args:
            target_atoms: Target atom group
            cutoff_distance: Cutoff distance in Angstroms
            
        Returns:
            Set of (chain_id, resid) tuples within cutoff
        """
        nearby_resids = set()
        
        # Get all protein residues (excluding target if it's protein)
        protein_residues = self.universe.select_atoms("protein").residues
        
        # Use cutoff directly since we specify it in Angstroms
        # Protein coordinates are typically in Angstroms
        cutoff_units = cutoff_distance  # cutoff_distance is already in Å
        print(f"Using improved distance calculation with cutoff: {cutoff_distance} Å")
        
        # Calculate distances from target atoms to residue atoms
        for residue in protein_residues:
            # Skip if same residue as target
            if (hasattr(target_atoms, 'residues') and 
                any(target_res.resid == residue.resid for target_res in target_atoms.residues)):
                continue
            
            # Check if any atom in this residue is within cutoff of any target atom
            residue_has_contact = False
            
            for target_atom in target_atoms:
                target_pos = target_atom.position
                
                for residue_atom in residue.atoms:
                    distance = np.linalg.norm(target_pos - residue_atom.position)
                    
                    if distance <= cutoff_units:
                        residue_has_contact = True
                        break
                
                if residue_has_contact:
                    break
            
            # If this residue has any contact, add it to nearby residues
            if residue_has_contact:
                # Get chain information
                chain_id = getattr(residue.atoms[0], 'chainID', 'Unknown')
                resid = residue.resid
                nearby_resids.add((chain_id, resid))
        
        return nearby_resids
    
    def find_nearby_residues_trajectory(self, target_atoms: mda.AtomGroup,
                                      cutoff_distance: float, 
                                      occupancy_threshold: float) -> Set[tuple]:
        """
        Find nearby residues using trajectory analysis with improved distance calculation
        
        Args:
            target_atoms: Target atom group
            cutoff_distance: Cutoff distance in Angstroms
            occupancy_threshold: Minimum occupancy fraction (0-1)
            
        Returns:
            Set of (chain_id, resid) tuples meeting occupancy threshold
        """
        if not self.trajectory_file:
            raise StructureAnalysisError("Trajectory file required for dynamic analysis")
        
        # Dictionary to track residue contact counts
        residue_contacts = {}
        total_frames = 0
        
        protein_residues = self.universe.select_atoms("protein").residues
        
        # Use cutoff directly since we specify it in Angstroms
        # Protein coordinates are typically in Angstroms
        cutoff_units = cutoff_distance  # cutoff_distance is already in Å
        print(f"Using improved distance calculation with cutoff: {cutoff_distance} Å for trajectory analysis")
        
        # Analyze trajectory
        for ts in self.universe.trajectory:
            total_frames += 1
            frame_contacts = set()
            
            # Find contacts in this frame using improved distance calculation
            for residue in protein_residues:
                # Skip if same residue as target
                if (hasattr(target_atoms, 'residues') and 
                    any(target_res.resid == residue.resid for target_res in target_atoms.residues)):
                    continue
                
                # Check if any atom in this residue is within cutoff of any target atom
                residue_has_contact = False
                
                for target_atom in target_atoms:
                    target_pos = target_atom.position
                    
                    for residue_atom in residue.atoms:
                        distance = np.linalg.norm(target_pos - residue_atom.position)
                        
                        if distance <= cutoff_units:
                            residue_has_contact = True
                            break
                    
                    if residue_has_contact:
                        break
                
                # If this residue has any contact, add it to frame contacts
                if residue_has_contact:
                    # Get chain information
                    chain_id = getattr(residue.atoms[0], 'chainID', 'Unknown')
                    resid = residue.resid
                    frame_contacts.add((chain_id, resid))
            
            # Update contact counts
            for residue_key in frame_contacts:
                residue_contacts[residue_key] = residue_contacts.get(residue_key, 0) + 1
        
        # Filter by occupancy threshold
        nearby_resids = set()
        for residue_key, count in residue_contacts.items():
            occupancy = count / total_frames
            if occupancy >= occupancy_threshold:
                nearby_resids.add(residue_key)
        
        return nearby_resids
    
    def find_nearby_residues_by_centroid(self, target_atoms: mda.AtomGroup, 
                                       cutoff_distance: float) -> Set[tuple]:
        """
        Find nearby residues using residue centroid distances (more accurate)
        
        Args:
            target_atoms: Target atom group
            cutoff_distance: Cutoff distance in Angstroms
            
        Returns:
            Set of (chain_id, resid) tuples within cutoff
        """
        nearby_resids = set()
        
        # Get all protein residues (excluding target if it's protein)
        protein_residues = self.universe.select_atoms("protein").residues
        
        # Convert cutoff_distance to nm if needed
        cutoff_nm = cutoff_distance / 10.0
        print(f"Using centroid-based cutoff distance: {cutoff_distance} Å ({cutoff_nm} nm)")
        
        # Calculate target centroid
        target_centroid = np.mean([atom.position for atom in target_atoms], axis=0)
        
        # Calculate distances from target centroid to residue centroids
        for residue in protein_residues:
            # Skip if same residue as target
            if (hasattr(target_atoms, 'residues') and 
                any(target_res.resid == residue.resid for target_res in target_atoms.residues)):
                continue
            
            # Calculate residue centroid
            residue_centroid = np.mean([atom.position for atom in residue.atoms], axis=0)
            
            # Calculate distance between centroids
            distance = np.linalg.norm(target_centroid - residue_centroid)
            
            if distance <= cutoff_nm:
                # Get chain information
                chain_id = getattr(residue.atoms[0], 'chainID', 'Unknown')
                resid = residue.resid
                nearby_resids.add((chain_id, resid))
        
        return nearby_resids
    
    def analyze_target_and_environment(self, target_selection: str, 
                                     cutoff_distance: float,
                                     use_trajectory: bool = False,
                                     occupancy_threshold: float = 0.5) -> Dict:
        """
        Complete analysis of target region and nearby residues using improved distance calculation
        
        Args:
            target_selection: Target selection string
            cutoff_distance: Cutoff distance in Angstroms
            use_trajectory: Whether to use trajectory for analysis
            occupancy_threshold: Occupancy threshold for trajectory analysis
            
        Returns:
            Dictionary with analysis results
        """
        # Store analysis parameters for later use
        self.target_selection = target_selection
        self.cutoff_distance = cutoff_distance
        self.use_trajectory = use_trajectory
        self.occupancy_threshold = occupancy_threshold
        
        # Identify target atoms
        target_atoms = self.identify_target_region(target_selection)
        
        # Find nearby residues using improved distance calculation
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
            'analysis_method': 'trajectory' if use_trajectory else 'static',
            'distance_method': 'improved_atom_to_atom'  # Updated method description
        }
        
        if use_trajectory:
            results['total_frames'] = len(self.universe.trajectory)
            results['occupancy_threshold'] = occupancy_threshold
        
        return results
    
    def _create_solute_selection(self, target_atoms: mda.AtomGroup, 
                               nearby_resids: Set[tuple]) -> mda.AtomGroup:
        """
        Create combined solute selection from target and nearby residues
        
        Args:
            target_atoms: Target atom group
            nearby_resids: Set of (chain_id, resid) tuples
            
        Returns:
            Combined atom group for solute
        """
        print(f"\nCreating solute selection:")
        print(f"  Target atom count: {len(target_atoms)}")
        print(f"  Nearby residue count: {len(nearby_resids)}")
        
        # Method 1: Try to create a combined selection string
        try:
            # Build selection parts
            selection_parts = []
            
            # Add target atoms (already have chain information from target_selection)
            if len(target_atoms) > 0:
                # Use the original target selection that was successful
                target_chain_id = getattr(target_atoms[0], 'chainID', None)
                if target_chain_id and target_chain_id != 'Unknown':
                    print(f"  Target selection: chainid {target_chain_id}")
                    selection_parts.append(f"chainid {target_chain_id}")
                else:
                    print(f"  Target selection: using target atom indices")
                    selection_parts.append(f"index {' '.join(map(str, target_atoms.indices))}")
            
            # Add nearby residues
            if nearby_resids:
                nearby_selections = []
                print(f"  Nearby residue selection:")
                for chain_id, resid in nearby_resids:
                    if chain_id != 'Unknown':
                        selection = f"chainid {chain_id} and resid {resid}"
                        nearby_selections.append(selection)
                        print(f"    {selection}")
                    else:
                        selection = f"resid {resid}"
                        nearby_selections.append(selection)
                        print(f"    {selection} (no chain info)")
                
                if nearby_selections:
                    selection_parts.append(f"({' or '.join(nearby_selections)})")
            
            # Combine selections
            if len(selection_parts) > 1:
                combined_selection = " or ".join([f"({part})" for part in selection_parts])
            else:
                combined_selection = selection_parts[0]
            
            print(f"  Combined selector: {combined_selection}")
            
            # Try the combined selection
            solute_atoms = self.universe.select_atoms(combined_selection)
            print(f"  ✓ Successfully created solute selection: {len(solute_atoms)} atoms")
            return solute_atoms
            
        except Exception as e:
            print(f"  ⚠ Combined selector failed: {e}")
            print(f"  Using fallback method...")
        
        # Method 2: Fallback - combine atom groups directly
        try:
            print(f"  Fallback method: directly combining atom groups")
            
            # Start with target atoms
            result_atoms = target_atoms.copy()
            
            # Add nearby residue atoms
            for chain_id, resid in nearby_resids:
                try:
                    if chain_id != 'Unknown':
                        nearby_atoms = self.universe.select_atoms(f"chainid {chain_id} and resid {resid}")
                    else:
                        nearby_atoms = self.universe.select_atoms(f"resid {resid}")
                    
                    if len(nearby_atoms) > 0:
                        result_atoms = result_atoms + nearby_atoms
                        print(f"    ✓ Added {chain_id}:{resid} - {len(nearby_atoms)} atoms")
                    else:
                        print(f"    ⚠ Not found {chain_id}:{resid}")
                        
                except Exception as e:
                    print(f"    ⚠ Selection {chain_id}:{resid} failed: {e}")
                    continue
            
            print(f"  ✓ Fallback method successful: {len(result_atoms)} atoms")
            return result_atoms
            
        except Exception as e:
            print(f"  ⚠ Fallback method also failed: {e}")
            print(f"  ✓ Using only target atoms: {len(target_atoms)} atoms")
            return target_atoms
    
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
        
        for residue_key in sorted(nearby_residues):
            try:
                # Handle both old format (int) and new format (tuple)
                if isinstance(residue_key, tuple):
                    chain_id, resid = residue_key
                    # Use chain-specific selection
                    if chain_id != 'Unknown':
                        residue = self.universe.select_atoms(f"chainid {chain_id} and resid {resid}").residues[0]
                    else:
                        residue = self.universe.select_atoms(f"resid {resid}").residues[0]
                else:
                    # Old format - just resid
                    resid = residue_key
                    chain_id = 'Unknown'
                    residue = self.universe.select_atoms(f"resid {resid}").residues[0]
                
                # Try different ways to get chain information
                actual_chain_id = chain_id
                # First try to get chain info from atoms (more reliable)
                if hasattr(residue.atoms[0], 'chainID') and residue.atoms[0].chainID:
                    actual_chain_id = residue.atoms[0].chainID
                elif hasattr(residue.atoms[0], 'segid') and residue.atoms[0].segid:
                    actual_chain_id = residue.atoms[0].segid
                # Then try residue level
                elif hasattr(residue, 'chainID') and residue.chainID:
                    actual_chain_id = residue.chainID
                elif hasattr(residue, 'segid') and residue.segid:
                    actual_chain_id = residue.segid
                elif hasattr(residue, 'chain') and residue.chain:
                    actual_chain_id = residue.chain
                
                info = {
                    'resid': resid,
                    'resname': residue.resname,
                    'chain': actual_chain_id,
                    'atom_count': len(residue.atoms)
                }
                residue_info.append(info)
            except (IndexError, AttributeError) as e:
                # Handle cases where residue might not be found
                if isinstance(residue_key, tuple):
                    chain_id, resid = residue_key
                else:
                    resid = residue_key
                    chain_id = 'Unknown'
                
                residue_info.append({
                    'resid': resid,
                    'resname': 'Unknown',
                    'chain': chain_id,
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
        
        # Convert nearby_residues to the format expected by solute_selector
        # IMPORTANT: Preserve chain information to avoid marking wrong chains
        nearby_residue_data = []
        for residue_key in nearby_residues:
            if isinstance(residue_key, tuple):
                chain_id, resid = residue_key
                nearby_residue_data.append({
                    'chain_id': chain_id,
                    'resid': resid
                })
            else:
                # Old format - assume unknown chain
                nearby_residue_data.append({
                    'chain_id': 'Unknown',
                    'resid': residue_key
                })
        
        return {
            'target_atom_indices': target_atoms.indices,
            'nearby_residue_data': nearby_residue_data,  # Now includes chain information
            'solute_atom_indices': solute_atoms.indices,
            'universe': self.universe
        }
    
    def print_selected_residues(self, analysis_results: Optional[Dict] = None, 
                               output_file: Optional[str] = None) -> None:
        """
        Print selected residues for user review and optionally save to file
        
        Args:
            analysis_results: Analysis results (uses stored results if None)
            output_file: Optional file path to save the residue information
        """
        residue_info = self.get_selected_residues_info(analysis_results)
        
        # Prepare the output content
        output_lines = []
        output_lines.append("="*60)
        output_lines.append("Selected Residues for REST2 Enhanced Sampling")
        output_lines.append("="*60)
        output_lines.append(f"{'ResID':<8} {'ResName':<8} {'Chain':<8} {'Atoms':<8}")
        output_lines.append("-"*60)
        
        for info in residue_info:
            line = f"{info['resid']:<8} {info['resname']:<8} {info['chain']:<8} {info['atom_count']:<8}"
            output_lines.append(line)
        
        output_lines.append("-"*60)
        output_lines.append(f"Total selected residues: {len(residue_info)}")
        
        if analysis_results:
            output_lines.append(f"Target atoms: {analysis_results['target_atom_count']}")
            output_lines.append(f"Total solute atoms: {analysis_results['solute_atom_count']}")
            output_lines.append(f"Analysis method: {analysis_results['analysis_method']}")
        
        output_lines.append("="*60)
        
        # Print to terminal
        for line in output_lines:
            print(line)
        print()  # Add extra newline for terminal output
        
        # Save to file if specified
        if output_file:
            try:
                output_path = Path(output_file)
                output_path.parent.mkdir(parents=True, exist_ok=True)
                
                with open(output_path, 'w') as f:
                    for line in output_lines:
                        f.write(line + '\n')
                
                print(f"✓ Residue information saved to: {output_path}")
                
                # Also save detailed information in a more structured format
                detailed_file = output_path.with_suffix('.detailed.txt')
                with open(detailed_file, 'w') as f:
                    f.write("Detailed Residue Information for REST2\n")
                    f.write("="*50 + "\n\n")
                    f.write(f"Analysis Summary:\n")
                    f.write(f"  Total selected residues: {len(residue_info)}\n")
                    if analysis_results:
                        f.write(f"  Target atoms: {analysis_results['target_atom_count']}\n")
                        f.write(f"  Total solute atoms: {analysis_results['solute_atom_count']}\n")
                        f.write(f"  Analysis method: {analysis_results['analysis_method']}\n")
                    f.write(f"  Cutoff distance: {getattr(self, 'cutoff_distance', 'N/A')}\n")
                    f.write(f"  Target selection: {getattr(self, 'target_selection', 'N/A')}\n\n")
                    
                    f.write("Residue Details:\n")
                    f.write("-"*50 + "\n")
                    f.write(f"{'ResID':<8} {'ResName':<8} {'Chain':<8} {'Atoms':<8} {'Distance':<10}\n")
                    f.write("-"*50 + "\n")
                    
                    # Add distance information if available
                    for info in residue_info:
                        distance_info = "N/A"
                        if hasattr(self, 'target_atoms') and self.target_atoms is not None:
                            try:
                                # Calculate distance from target to this residue
                                residue_atoms = self.universe.select_atoms(f"resid {info['resid']}")
                                if len(residue_atoms) > 0:
                                    min_distance = float('inf')
                                    for target_atom in self.target_atoms:
                                        for residue_atom in residue_atoms:
                                            dist = np.linalg.norm(target_atom.position - residue_atom.position)
                                            min_distance = min(min_distance, dist)
                                    if min_distance != float('inf'):
                                        distance_info = f"{min_distance:.2f}"
                            except:
                                distance_info = "N/A"
                        
                        f.write(f"{info['resid']:<8} {info['resname']:<8} {info['chain']:<8} "
                               f"{info['atom_count']:<8} {distance_info:<10}\n")
                    
                    f.write("-"*50 + "\n")
                    f.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                
                print(f"✓ Detailed information saved to: {detailed_file}")
                
            except Exception as e:
                print(f"Warning: Could not save to file {output_file}: {e}")

    def debug_distance_calculation(self, target_selection: str, cutoff_distance: float,
                                 max_debug_residues: int = 10) -> None:
        """
        Debug method to verify distance calculations
        
        Args:
            target_selection: Target selection string
            cutoff_distance: Cutoff distance in Angstroms
            max_debug_residues: Maximum number of residues to show in debug output
        """
        print(f"\n=== Distance Calculation Debug ===")
        print(f"Target selection: {target_selection}")
        print(f"Cutoff distance: {cutoff_distance} Å")
        
        # Identify target atoms
        target_atoms = self.identify_target_region(target_selection)
        print(f"Target atoms found: {len(target_atoms)}")
        
        # Convert cutoff to nm
        cutoff_nm = cutoff_distance / 10.0
        if cutoff_distance > 10:
            cutoff_nm = cutoff_distance / 10.0
            print(f"Converted cutoff: {cutoff_nm} nm")
        
        # Get protein residues
        protein_residues = self.universe.select_atoms("protein").residues
        print(f"Total protein residues: {len(protein_residues)}")
        
        # Calculate target centroid
        target_centroid = np.mean([atom.position for atom in target_atoms], axis=0)
        print(f"Target centroid: {target_centroid}")
        
        # Show distances to first few residues
        print(f"\nDistances to first {max_debug_residues} residues:")
        print(f"{'ResID':<8} {'ResName':<8} {'Chain':<8} {'Distance (nm)':<12} {'Within cutoff':<12}")
        print("-" * 60)
        
        count = 0
        for residue in protein_residues:
            if count >= max_debug_residues:
                break
                
            # Skip if same residue as target
            if (hasattr(target_atoms, 'residues') and 
                any(target_res.resid == residue.resid for target_res in target_atoms.residues)):
                continue
            
            # Calculate residue centroid
            residue_centroid = np.mean([atom.position for atom in residue.atoms], axis=0)
            distance = np.linalg.norm(target_centroid - residue_centroid)
            
            # Get chain info
            chain_id = getattr(residue.atoms[0], 'chainID', 'Unknown')
            
            within_cutoff = "✓" if distance <= cutoff_nm else "✗"
            
            print(f"{residue.resid:<8} {residue.resname:<8} {chain_id:<8} {distance:<12.3f} {within_cutoff:<12}")
            count += 1
        
        print("-" * 60)
        print("=== End Debug ===\n")


def main():
    """Test structure analyzer"""
    try:
        # Example usage
        structure_file = "example/MD_results/md.gro"
        topology_file = "example/MD_results/md.tpr"  # TPR file for chain information
        trajectory_file = "example/MD_results/md.xtc"  # Optional
        
        # Initialize analyzer
        analyzer = StructureAnalyzer(structure_file, topology_file, trajectory_file)
        
        # Debug distance calculation first
        print("=== Testing Distance Calculation ===")
        analyzer.debug_distance_calculation("chainid C", 4.0)
        
        # Analyze peptide target with centroid-based distance
        print("=== Analysis with Centroid-based Distance ===")
        results = analyzer.analyze_target_and_environment(
            target_selection="chainid C",
            cutoff_distance=4.0,
            use_trajectory=False
        )
        
        # Print selected residues for user review
        analyzer.print_selected_residues(results)
        
        # Get data for next module
        solute_data = analyzer.get_solute_selection_data(results)
        print(f"\nData prepared for solute selector:")
        print(f"  - Target atoms: {len(solute_data['target_atom_indices'])}")
        print(f"  - Nearby residues: {len(solute_data['nearby_residue_data'])}")
        print(f"  - Total solute atoms: {len(solute_data['solute_atom_indices'])}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()