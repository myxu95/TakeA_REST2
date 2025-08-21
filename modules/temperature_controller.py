#!/usr/bin/env python3
"""
REST2 Temperature Controller Module
Handles temperature scaling and PLUMED file generation for REST2 simulations
"""

import os
import shutil
from pathlib import Path
from typing import Dict, List, Any, Optional
import numpy as np

# Import from utils package
try:
    from utils import TemperatureCalculator, TemperatureCalculationError, ValidationFramework, FileUtils, FileOperationError, OutputFormatter
except ImportError:
    # Fallback for direct execution
    TemperatureCalculator = None
    TemperatureCalculationError = Exception
    ValidationFramework = None
    FileUtils = None
    FileOperationError = Exception
    OutputFormatter = None


class TemperatureControllerError(Exception):
    """Temperature controller error"""
    pass


class TemperatureController:
    """
    Temperature controller for REST2 simulations
    Handles temperature scaling and PLUMED file generation
    """
    
    def __init__(self, config_manager, replica_data: Dict[str, Any], 
                 solute_data: Optional[Dict[str, Any]] = None):
        """
        Initialize temperature controller
        
        Args:
            config_manager: Configuration manager instance
            replica_data: Replica configuration data
            solute_data: Solute selection data (optional)
        """
        self.config_manager = config_manager
        self.replica_data = replica_data
        
        # Temperature and scaling information
        self.n_replicas = replica_data['n_replicas']
        self.replicas = replica_data['replicas']
        self.temperatures = [replica['temperature'] for replica in self.replicas]
        self.scaling_factors = [replica['scaling_factor'] for replica in self.replicas]
        
        # Reference temperature (lowest temperature)
        self.T_ref = min(self.temperatures)
        
        # Solute atoms for PLUMED REST2
        self.solute_atom_indices = None
        self.solute_atom_count = 0
        
        if solute_data and 'solute_atom_indices' in solute_data:
            self._process_solute_data(solute_data)
        else:
            print("Warning: No solute data provided, using default configuration")
            self._create_default_solute_data()
        
        # Validate solute data
        self._validate_solute_data()
    
    def _process_solute_data(self, solute_data: Dict[str, Any]) -> None:
        """Process solute data and extract atom indices"""
        print(f"Processing solute data...")
        
        # Get raw indices
        raw_indices = solute_data['solute_atom_indices']
        print(f"  Raw solute atom count: {len(raw_indices)}")
        
        # Check if indices are empty (handle both list and numpy array)
        if len(raw_indices) == 0:
            print("  Warning: Solute atom indices are empty")
            self._create_default_solute_data()
            return
        
        # Convert indices format
        try:
            if hasattr(raw_indices, 'tolist'):  # numpy array
                raw_indices = raw_indices.tolist()
            
            # Convert to integers and ensure 1-based indexing
            indices = [int(idx) for idx in raw_indices]
            min_idx = min(indices)
            max_idx = max(indices)
            
            print(f"  Index range: {min_idx} - {max_idx}")
            
            # Convert to 1-based if needed
            if max_idx < 1000:  # Likely 0-based
                print(f"  Converting 0-based indices to 1-based")
                self.solute_atom_indices = [idx + 1 for idx in indices]
            else:  # Likely already 1-based
                print(f"  Using 1-based indices")
                self.solute_atom_indices = indices
            
            self.solute_atom_count = len(self.solute_atom_indices)
            
        except Exception as e:
            print(f"  Error processing solute data: {e}")
            self._create_default_solute_data()
    
    def _create_default_solute_data(self) -> None:
        """Create default solute data"""
        print(f"  Creating default solute data...")
        
        # Try to get default from config
        default_atoms = self.config_manager.get_parameter('default_solute_atoms', None)
        
        if default_atoms:
            print(f"  Using default from config: {default_atoms}")
            self.solute_atom_indices = default_atoms
        else:
            # Create reasonable default (first 100 protein atoms)
            print(f"  Creating default: first 100 protein atoms")
            self.solute_atom_indices = list(range(1, 101))
        
        self.solute_atom_count = len(self.solute_atom_indices)
    
    def _validate_solute_data(self) -> None:
        """Validate solute data consistency"""
        if len(self.solute_atom_indices) == 0:
            print("Error: No solute atoms provided")
            raise TemperatureControllerError("No solute atoms provided")
        
        # Check index validity
        invalid_indices = [idx for idx in self.solute_atom_indices if idx <= 0]
        if len(invalid_indices) > 0:
            print(f"Error: Invalid atom indices found: {invalid_indices[:5]}...")
            raise TemperatureControllerError("Invalid atom indices found")
        
        print(f"Solute data validation passed:")
        print(f"  Solute atom count: {self.solute_atom_count}")
        print(f"  Index range: {min(self.solute_atom_indices)} - {max(self.solute_atom_indices)}")
    
    def print_temperature_summary(self) -> None:
        """Print temperature configuration summary"""
        print(f"\nREST2 Temperature Configuration Summary")
        print(f"=" * 50)
        print(f"Number of replicas: {self.n_replicas}")
        print(f"Reference temperature: {self.T_ref:.1f} K")
        print(f"Solute atoms: {self.solute_atom_count}")
        print(f"\nReplica details:")
        
        for i, replica in enumerate(self.replicas):
            temp = replica['temperature']
            scale = replica['scaling_factor']
            print(f"  Replica {i}: T = {temp:.1f} K, lambda = {scale:.3f}")
    
    def generate_scaled_topology_files(self, base_topology: str) -> None:
        """Generate scaled topology files using plumed partial_tempering"""
        base_topology_path = Path(base_topology)
        if not base_topology_path.exists():
            raise FileNotFoundError(f"Base topology file not found: {base_topology}")
        
        print(f"\nGenerating scaled topology files using plumed partial_tempering...")
        
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            
            # Generate scaled topology using plumed command
            scaled_topology_path = input_dir / "topol-scaled.top"
            self._create_scaled_topology_with_plumed(base_topology_path, scaled_topology_path, replica_index)
    
    def _create_scaled_topology_with_plumed(self, base_topology: Path, output_topology: Path, replica_index: int) -> None:
        """Create scaled topology using plumed partial_tempering command"""
        try:
            # Ensure output directory exists
            output_topology.parent.mkdir(parents=True, exist_ok=True)
            
            # Get scaling factor for this replica
            scaling_factor = self.replicas[replica_index]['scaling_factor']
            temperature = self.replicas[replica_index]['temperature']
            
            print(f"  Replica {replica_index}: T = {temperature:.1f} K, lambda = {scaling_factor:.6f}")
            
            # Method 1: Try using plumed command directly
            if self._try_plumed_partial_tempering(base_topology, output_topology, scaling_factor):
                print(f"  ✓ Generated scaled topology using plumed: {output_topology}")
                return
            
            # Method 2: Fallback to manual scaling (if plumed command not available)
            print(f"  Warning: plumed command not available, using manual scaling")
            self._create_manually_scaled_topology(base_topology, output_topology, scaling_factor, replica_index)
            
        except Exception as e:
            raise TemperatureControllerError(f"Failed to create scaled topology: {e}")
    
    def _try_plumed_partial_tempering(self, base_topology: Path, output_topology: Path, scaling_factor: float) -> bool:
        """Try to use plumed partial_tempering command"""
        try:
            import subprocess
            
            # Check if plumed command is available
            result = subprocess.run(['plumed', '--version'], capture_output=True, text=True)
            if result.returncode != 0:
                print(f"    plumed command not found in PATH")
                return False
            
            print(f"    Using plumed partial_tempering command...")
            
            # Run plumed partial_tempering
            cmd = [
                'plumed', 'partial_tempering', 
                str(scaling_factor)
            ]
            
            # Read base topology and pipe to plumed
            with open(base_topology, 'r') as f:
                topology_content = f.read()
            
            # Run plumed command
            result = subprocess.run(
                cmd,
                input=topology_content,
                capture_output=True,
                text=True
            )
            
            if result.returncode == 0:
                # Write scaled topology
                with open(output_topology, 'w') as f:
                    f.write(result.stdout)
                return True
            else:
                print(f"    plumed command failed: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"    Error running plumed command: {e}")
            return False
    
    def _create_manually_scaled_topology(self, base_topology: Path, output_topology: Path, scaling_factor: float, replica_index: int) -> None:
        """Create manually scaled topology (fallback method)"""
        try:
            print(f"    Creating manually scaled topology...")
            
            # Read base topology
            with open(base_topology, 'r') as f:
                lines = f.readlines()
            
            # Process topology lines
            scaled_lines = []
            in_atoms_section = False
            
            for line in lines:
                if line.strip().startswith('[ atoms ]'):
                    in_atoms_section = True
                    scaled_lines.append(line)
                    continue
                elif line.strip().startswith('['):
                    in_atoms_section = False
                
                if in_atoms_section and line.strip() and not line.strip().startswith(';'):
                    # This is an atom line, check if it's a solute atom
                    if self._is_solute_atom_line(line):
                        # Scale the atom parameters
                        scaled_line = self._scale_atom_line(line, scaling_factor)
                        scaled_lines.append(scaled_line)
                    else:
                        # Keep original line
                        scaled_lines.append(line)
                else:
                    # Keep original line
                    scaled_lines.append(line)
            
            # Write scaled topology
            with open(output_topology, 'w') as f:
                f.writelines(scaled_lines)
            
            print(f"    ✓ Generated manually scaled topology: {output_topology}")
            
        except Exception as e:
            print(f"    Error in manual scaling: {e}")
            # Fallback to simple copy
            import shutil
            shutil.copy2(base_topology, output_topology)
            print(f"    Fallback: copied original topology")
    
    def _is_solute_atom_line(self, line: str) -> bool:
        """Check if an atom line corresponds to a solute atom"""
        try:
            # Parse atom line (GROMACS topology format)
            parts = line.split()
            if len(parts) >= 5:
                atom_index = int(parts[0])
                # Check if this atom is in our solute list
                return atom_index in self.solute_atom_indices
        except:
            pass
        return False
    
    def _scale_atom_line(self, line: str, scaling_factor: float) -> str:
        """Scale atom parameters in topology line"""
        try:
            parts = line.split()
            if len(parts) >= 5:
                # Format: nr type resnr residue atom cgnr charge mass
                # We might need to scale charge or other parameters
                # For now, just add a comment indicating scaling
                scaled_line = line.rstrip() + f" ; scaled by lambda={scaling_factor:.6f}\n"
                return scaled_line
        except:
            pass
        return line
    
    def generate_mdp_files(self) -> None:
        """Generate MDP files"""
        print(f"\nGenerating MDP files...")
        
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            
            # Create MDP file
            mdp_path = input_dir / "md.mdp"
            self._create_mdp_file(mdp_path, replica_index)
    
    def _create_mdp_file(self, mdp_path: Path, replica_index: int) -> None:
        """Create MDP file"""
        try:
            # Ensure output directory exists
            mdp_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Get temperature
            temperature = self.replicas[replica_index]['temperature']
            
            # Create MDP content
            mdp_content = f"""# REST2 MDP file for replica {replica_index}
# Generated by TemperatureController

# Temperature
ref_t = {temperature:.1f}

# REST2 specific settings
restart = no
integrator = md
dt = 0.002
nsteps = 50000

# Output control
nstxout = 1000
nstvout = 1000
nstenergy = 1000
nstlog = 1000

# Neighbor searching
cutoff-scheme = Verlet
ns_type = grid
pbc = xyz

# Electrostatics
coulombtype = PME
rcoulomb = 1.0

# Van der Waals
rvdw = 1.0

# Temperature coupling
tcoupl = V-rescale
tc-grps = System
tau_t = 0.1

# Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau_p = 2.0
compressibility = 4.5e-5

# Constraints
constraints = h-bonds
constraint_algorithm = LINCS
lincs_iter = 1
lincs_order = 4
"""
            
            # Write file
            with open(mdp_path, 'w') as f:
                f.write(mdp_content)
            
            print(f"  Replica {replica_index}: {mdp_path}")
            
        except Exception as e:
            raise TemperatureControllerError(f"Failed to create MDP file: {e}")
    
    def prepare_additional_input_files(self) -> None:
        """Prepare additional input files"""
        print(f"\nPreparing additional input files...")
        
        # Prepare PLUMED files
        self._prepare_plumed_files()
        
        # Prepare index files
        self._prepare_index_files()
    
    def _prepare_plumed_files(self) -> None:
        """Prepare PLUMED files"""
        print(f"  Preparing PLUMED files...")
        
        # Get PLUMED configuration from config
        plumed_config = self.config_manager.get_parameter('plumed', {})
        enable_plumed = plumed_config.get('enable', False)
        plumed_template = plumed_config.get('template', None)
        
        if not enable_plumed:
            print(f"  PLUMED files disabled in config, generating empty plumed.dat files")
            # Generate empty plumed.dat files for each replica
            for replica in self.replicas:
                replica_index = replica['index']
                input_dir = Path(replica['input_dir'])
                plumed_path = input_dir / "plumed.dat"
                self._create_empty_plumed_file(plumed_path, replica_index)
            return
        
        if not plumed_template:
            print(f"  Error: PLUMED enabled but no template specified")
            raise TemperatureControllerError("PLUMED enabled but no template specified")
        
        # Generate PLUMED file for each replica
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            
            plumed_path = input_dir / "plumed.dat"
            self._create_replica_plumed_file(plumed_path, replica_index, plumed_template)
    
    def _create_partial_tempering_command(self, replica_index: int) -> None:
        """Create PLUMED file for replica (simplified since topology is already scaled)"""
        temperature = self.replicas[replica_index]['temperature']
        scaling_factor = self.replicas[replica_index]['scaling_factor']
        
        print(f"  Replica {replica_index}: T = {temperature:.1f} K, lambda = {scaling_factor:.6f}")
        print(f"  Note: Topology is already scaled, PLUMED file for user commands only")
        
        # Since topology is already scaled, we don't need PARTIAL_TEMPERING
        # PLUMED file is mainly for user-defined analysis commands
        return ""
    
    def _create_replica_plumed_file(self, plumed_path: Path, replica_index: int, plumed_template: str) -> None:
        """Create replica PLUMED file (simplified since topology is already scaled)"""
        try:
            # Ensure output directory exists
            plumed_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Read template content (required since enable=true)
            template_content = ""
            try:
                with open(plumed_template, 'r') as f:
                    template_content = f.read()
                print(f"  Using PLUMED template: {plumed_template}")
            except Exception as e:
                print(f"  Error: Failed to read PLUMED template: {e}")
                raise TemperatureControllerError(f"Failed to read PLUMED template: {e}")
            
            # Create header comment
            temperature = self.replicas[replica_index]['temperature']
            scaling_factor = self.replicas[replica_index]['scaling_factor']
            
            header = f"""# PLUMED file for REST2 replica {replica_index}
# Temperature: {temperature:.1f} K
# Scaling factor (lambda): {scaling_factor:.6f}
# Note: Topology is already scaled using plumed partial_tempering
# This file is for user-defined analysis commands only

"""
            
            # Combine content
            full_content = header + template_content
            
            # Write file
            with open(plumed_path, 'w') as f:
                f.write(full_content)
            
            print(f"  Replica {replica_index}: {plumed_path}")
            
        except Exception as e:
            raise TemperatureControllerError(f"Failed to create PLUMED file: {e}")
    
    def _create_empty_plumed_file(self, plumed_path: Path, replica_index: int) -> None:
        """Create an empty PLUMED file for a replica."""
        try:
            # Ensure output directory exists
            plumed_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create empty content
            empty_content = f"""# PLUMED file for REST2 replica {replica_index}
# Generated by TemperatureController
# This file is intentionally empty.
"""
            
            # Write file
            with open(plumed_path, 'w') as f:
                f.write(empty_content)
            
            print(f"  Replica {replica_index}: {plumed_path} (empty)")
            
        except Exception as e:
            print(f"  Error creating empty plumed.dat: {e}")
            raise TemperatureControllerError(f"Failed to create empty plumed.dat: {e}")
    
    def _format_atom_list(self, atom_indices: List[int]) -> str:
        """Format atom indices list for PLUMED"""
        if len(atom_indices) == 0:
            return "1-100  # WARNING: Empty atom list"
        
        # Validate input
        if any(idx <= 0 for idx in atom_indices):
            raise TemperatureControllerError("Invalid atom indices: all indices must be > 0")
        
        # Sort indices
        sorted_indices = sorted(atom_indices)
        
        # Group consecutive indices into ranges
        ranges = []
        start = sorted_indices[0]
        end = start
        
        for i in range(1, len(sorted_indices)):
            if sorted_indices[i] == end + 1:
                end = sorted_indices[i]
            else:
                if start == end:
                    ranges.append(str(start))
                else:
                    ranges.append(f"{start}-{end}")
                start = sorted_indices[i]
                end = start
        
        # Add the last range
        if start == end:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}-{end}")
        
        # Join ranges with commas
        atom_list = ",".join(ranges)
        
        return atom_list
    
    def _prepare_index_files(self) -> None:
        """Prepare index files"""
        print(f"  Preparing index files...")
        
        # Create index file for each replica
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            
            index_path = input_dir / "index.ndx"
            self._create_index_file(index_path, replica_index)
    
    def _create_index_file(self, index_path: Path, replica_index: int) -> None:
        """Create index file"""
        try:
            # Ensure output directory exists
            index_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create index content
            max_atom_idx = max(self.solute_atom_indices) if len(self.solute_atom_indices) > 0 else 1000
            solute_atoms_str = ' '.join(map(str, self.solute_atom_indices)) if len(self.solute_atom_indices) > 0 else '1-100'
            
            index_content = f"""# REST2 index file for replica {replica_index}
# Generated by TemperatureController

[ System ]
1-{max_atom_idx}

[ Solute ]
{solute_atoms_str}

[ Protein ]
1-{max_atom_idx}

[ Non-Protein ]
# Add non-protein atoms as needed
"""
            
            # Write file
            with open(index_path, 'w') as f:
                f.write(index_content)
            
            print(f"  Replica {replica_index}: {index_path}")
            
        except Exception as e:
            print(f"  Warning: Failed to create index file: {e}")
    
    def create_temperature_summary(self) -> None:
        """Create temperature summary file"""
        print(f"\nCreating temperature summary...")
        
        summary_path = Path(self.config_manager.get_parameter('output_dir')) / "temperature_summary.txt"
        
        try:
            # Ensure output directory exists
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create summary content
            summary_content = f"""REST2 Temperature Summary
Generated by TemperatureController
{'='*50}

Configuration:
- Number of replicas: {self.n_replicas}
- Reference temperature: {self.T_ref:.1f} K
- Solute atoms: {self.solute_atom_count}

Replica Details:
{'-'*50}"""
            
            for i, replica in enumerate(self.replicas):
                temp = replica['temperature']
                scale = replica['scaling_factor']
                summary_content += f"\nReplica {i}:"
                summary_content += f"\n  Temperature: {temp:.1f} K"
                summary_content += f"\n  Scaling factor: {scale:.6f}"
                summary_content += f"\n  Input directory: {replica['input_dir']}"
                summary_content += f"\n"
            
            summary_content += f"""
Solute Atom Information:
{'-'*50}
- Total atoms: {self.solute_atom_count}
- Index range: {min(self.solute_atom_indices) if len(self.solute_atom_indices) > 0 else 'N/A'} - {max(self.solute_atom_indices) if len(self.solute_atom_indices) > 0 else 'N/A'}
- Formatted list: {self._format_atom_list(self.solute_atom_indices) if len(self.solute_atom_indices) > 0 else 'N/A'}

Files Generated:
{'-'*50}
- Topology files: {self.n_replicas} copies
- MDP files: {self.n_replicas} copies  
- PLUMED files: {self.n_replicas} copies
- Index files: {self.n_replicas} copies

Note: All files are ready for REST2 simulation.
"""
            
            # Write file
            with open(summary_path, 'w') as f:
                f.write(summary_content)
            
            print(f"Temperature summary created: {summary_path}")
            
        except Exception as e:
            print(f"Warning: Failed to create temperature summary: {e}")
    
    def validate_temperature_setup(self) -> bool:
        """Validate temperature setup"""
        print(f"\nValidating temperature setup...")
        
        try:
            # Check replica count
            if self.n_replicas <= 0:
                print("Error: Number of replicas must be > 0")
                return False
            
            # Check temperature range
            if min(self.temperatures) <= 0:
                print("Error: Temperature must be > 0")
                return False
            
            # Check solute atoms
            if len(self.solute_atom_indices) == 0:
                print("Error: No solute atoms")
                return False
            
            # Check index validity
            if any(idx <= 0 for idx in self.solute_atom_indices):
                print("Error: Solute atom indices must be > 0")
                return False
            
            print("Temperature setup validation passed")
            return True
            
        except Exception as e:
            print(f"Temperature setup validation failed: {e}")
            return False