#!/usr/bin/env python3
"""
REST2 Temperature Controller Module
Handles temperature-specific topology modifications and input file preparation
"""

import os
import shutil
import re
from pathlib import Path
from typing import Dict, List, Any, Optional
import numpy as np


class TemperatureControllerError(Exception):
    """Temperature controller error"""
    pass


class TemperatureController:
    """
    Temperature controller for REST2 simulations
    Generates temperature-specific topology files and input files
    """
    
    def __init__(self, config_manager, replica_data: Dict[str, Any], 
                 solute_data: Optional[Dict[str, Any]] = None):
        """
        Initialize temperature controller
        
        Args:
            config_manager: ConfigManager instance
            replica_data: Data from ReplicaGenerator containing replica information
            solute_data: Data from SoluteSelector containing solute atom information
        """
        self.config = config_manager
        self.replica_data = replica_data
        self.solute_data = solute_data
        self.replicas = replica_data['replicas']
        self.n_replicas = replica_data['n_replicas']
        
        # Temperature and scaling information
        self.temperatures = [replica['temperature'] for replica in self.replicas]
        self.scaling_factors = [replica['scaling_factor'] for replica in self.replicas]
        
        # Reference temperature (lowest temperature)
        self.T_ref = min(self.temperatures)
        
        # Solute atoms for PLUMED REST2
        self.solute_atom_indices = None
        if solute_data and 'solute_atom_indices' in solute_data:
            # Convert to 1-based indexing for PLUMED and ensure it's a Python list
            solute_indices = solute_data['solute_atom_indices']
            if hasattr(solute_indices, 'tolist'):  # numpy array
                self.solute_atom_indices = [int(idx) + 1 for idx in solute_indices.tolist()]
            else:  # Python list
                self.solute_atom_indices = [int(idx) + 1 for idx in solute_indices]
            
            # Validate solute data
            self._validate_solute_data(solute_data)
    
    def _validate_solute_data(self, solute_data: Dict[str, Any]) -> None:
        """Validate solute data consistency"""
        solute_atoms = solute_data.get('solute_atom_indices')
        if solute_atoms is None or len(solute_atoms) == 0:
            print("Warning: No solute atoms provided for REST2 scaling")
            return
        
        solute_count = len(solute_atoms)
        print(f"Solute atoms: {solute_count}")
        
        # Validate reasonable counts
        if solute_count < 10:
            print("Warning: Very few solute atoms (< 10)")
        elif solute_count > 10000:
            print("Warning: Very many solute atoms (> 10000)")
        
        # Check target atoms
        target_atoms = solute_data.get('target_atom_indices', [])
        if target_atoms is not None and len(target_atoms) > 0:
            target_count = len(target_atoms)
            print(f"Target atoms: {target_count}")
        
        # Check nearby residues
        nearby_residues = solute_data.get('nearby_residue_ids', [])
        if nearby_residues is not None and len(nearby_residues) > 0:
            nearby_count = len(nearby_residues)
            print(f"Nearby residues: {nearby_count}")
    
    def generate_scaled_topology_files(self, base_topology: str) -> None:
        """
        Generate topology files for each replica (no scaling needed for PLUMED REST2)
        
        Args:
            base_topology: Path to base REST2-modified topology file
        """
        base_topology_path = Path(base_topology)
        if not base_topology_path.exists():
            raise FileNotFoundError(f"Base topology file not found: {base_topology}")
        
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            
            # For PLUMED REST2, we just copy the base topology
            # The scaling is handled by PLUMED PARTIAL_TEMPERING
            scaled_topology_path = input_dir / "topol.top"
            self._create_replica_topology(
                base_topology_path, 
                scaled_topology_path, 
                replica_index
            )
    
    def _create_replica_topology(self, base_topology: Path, output_topology: Path,
                               replica_index: int) -> None:
        """
        Create replica topology file (for PLUMED REST2, no parameter scaling needed)
        
        Args:
            base_topology: Base topology file path
            output_topology: Output topology file path
            replica_index: Replica index
        """
        with open(base_topology, 'r') as f:
            content = f.read()
        
        # Write topology
        with open(output_topology, 'w') as f:
            f.write(content)
    
    def generate_mdp_files(self) -> None:
        """Generate MDP files for each replica with appropriate temperature settings"""
        base_mdp_content = self._create_base_mdp_template()
        
        for replica in self.replicas:
            replica_index = replica['index']
            temperature = replica['temperature']
            input_dir = Path(replica['input_dir'])
            
            # Customize MDP for this replica
            mdp_content = self._customize_mdp_for_replica(
                base_mdp_content, 
                replica_index, 
                temperature
            )
            
            # Write MDP file
            mdp_file = input_dir / "rest2.mdp"
            with open(mdp_file, 'w') as f:
                f.write(mdp_content)
    
    def _create_base_mdp_template(self) -> str:
        """Create base MDP template for REST2 simulation"""
        replex = self.config.get_parameter('replex')
        
        mdp_template = f"""title                   = REST2 Enhanced Sampling
; Run parameters
integrator              = md            ; leap-frog integrator
nsteps                  = NSTEPS_PLACEHOLDER    ; Will be set from config
dt                      = 0.002         ; 2 fs

; Output control
nstxout                 = 0             ; suppress bulky .trr file
nstvout                 = 0             ; suppress velocities output
nstfout                 = 0             ; suppress forces output
nstenergy               = 5000          ; save energies every 10.0 ps
nstlog                  = 5000          ; update log file every 10.0 ps
nstxout-compressed      = 5000          ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System        ; save the whole system

; Bond parameters
continuation            = yes           ; Restarting after equilibration
constraint_algorithm    = lincs         ; holonomic constraints
constraints             = h-bonds       ; bonds involving H are constrained
lincs_iter              = 1             ; accuracy of LINCS
lincs_order             = 4             ; also related to accuracy

; Neighbor searching
cutoff-scheme           = Verlet        ; Buffered neighbor searching
ns_type                 = grid          ; search neighboring grid cells
nstlist                 = 10            ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0           ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0           ; short-range van der Waals cutoff (in nm)

; Electrostatics
coulombtype             = PME           ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4             ; cubic interpolation
fourierspacing          = 0.16          ; grid spacing for FFT

; Temperature coupling
tcoupl                  = V-rescale     ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = TEMP_PLACEHOLDER TEMP_PLACEHOLDER   ; reference temperature

; Pressure coupling
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1

; Periodic boundary conditions
pbc                     = xyz           ; 3-D PBC

; Dispersion correction
DispCorr                = EnerPres      ; account for cut-off vdW scheme

; Velocity generation
gen_vel                 = no            ; Velocity generation is off

; Replica exchange
nstreplex               = {replex}
"""
        return mdp_template
    
    def _customize_mdp_for_replica(self, base_mdp: str, replica_index: int, 
                                 temperature: float) -> str:
        """
        Customize MDP content for specific replica
        
        Args:
            base_mdp: Base MDP template
            replica_index: Index of current replica
            temperature: Temperature for this replica
            
        Returns:
            Customized MDP content
        """
        # Calculate nsteps from simulation time
        production_time_ns = self.config.get_parameter('simulation.production_time', 100.0)
        dt_ps = 0.002  # 2 fs
        nsteps = int(production_time_ns * 1000 / dt_ps)  # Convert ns to steps
        
        # Replace placeholders
        customized_mdp = base_mdp.replace("TEMP_PLACEHOLDER", f"{temperature:.1f}")
        customized_mdp = customized_mdp.replace("NSTEPS_PLACEHOLDER", str(nsteps))
        
        # Add replica-specific header
        header = f"""; Replica {replica_index} MDP File
; Temperature: {temperature:.1f} K
; REST2 scaling factor (λ): {self.scaling_factors[replica_index]:.6f}
; Simulation time: {production_time_ns:.1f} ns ({nsteps} steps)

"""
        
        return header + customized_mdp
    
    def prepare_additional_input_files(self) -> None:
        """Prepare additional input files needed for each replica (不再生成index.ndx)"""
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            # 只处理PLUMED文件
            self._prepare_plumed_file(input_dir, replica_index)
    
    def _prepare_plumed_file(self, input_dir: Path, replica_index: int) -> None:
        """
        Prepare PLUMED file with REST2 PARTIAL_TEMPERING for replica
        
        Args:
            input_dir: Replica input directory
            replica_index: Replica index
        """
        plumed_dat = self.config.get_parameter('plumed_dat')
        
        # Add PARTIAL_TEMPERING command for REST2
        rest2_command = self._create_partial_tempering_command(replica_index)
        
        # Read original PLUMED content if exists
        original_content = ""
        if plumed_dat and Path(plumed_dat).exists():
            with open(plumed_dat, 'r') as f:
                original_content = f.read()
        
        # Write PLUMED file
        with open(input_dir / "plumed.dat", 'w') as f:
            f.write(rest2_command)
            f.write("\n")
            if original_content:
                f.write(original_content)
    
    def _create_partial_tempering_command(self, replica_index: int) -> str:
        """
        Create PLUMED PARTIAL_TEMPERING command for REST2
        
        Args:
            replica_index: Replica index
            
        Returns:
            PARTIAL_TEMPERING command string
        """
        temperature = self.temperatures[replica_index]
        scaling_factor = self.scaling_factors[replica_index]
        
        # Create atom list for PARTIAL_TEMPERING
        if self.solute_atom_indices is not None and len(self.solute_atom_indices) > 0:
            atom_list = self._format_atom_list(self.solute_atom_indices)
        else:
            atom_list = "1-100  # EDIT THIS: Replace with actual solute atom indices"
        
        partial_tempering_cmd = f"""PARTIAL_TEMPERING ...
  ATOMS={atom_list}
  TEMP={self.T_ref:.1f}
  LAMBDA={scaling_factor:.6f}
  LABEL=rest2_scaling
... PARTIAL_TEMPERING

"""
        
        return partial_tempering_cmd
    
    def _format_atom_list(self, atom_indices: List[int]) -> str:
        """
        Format atom indices list for PLUMED
        
        Args:
            atom_indices: List of atom indices (1-based)
            
        Returns:
            Formatted atom list string
        """
        if not atom_indices:
            return "1-100  # EDIT THIS: Replace with actual solute atom indices"
        
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
        
        # Join ranges with better formatting
        atom_list = ",".join(ranges)
        
        # If the list is too long, break it into multiple lines
        if len(atom_list) > 80:  # Increased threshold for better readability
            formatted_ranges = []
            current_line = ""
            line_length = 0
            
            for range_str in ranges:
                if line_length + len(range_str) + 1 > 80:  # +1 for comma
                    if current_line:
                        formatted_ranges.append(current_line.rstrip())
                    current_line = range_str
                    line_length = len(range_str)
                else:
                    if current_line:
                        current_line += "," + range_str
                        line_length += len(range_str) + 1
                    else:
                        current_line = range_str
                        line_length = len(range_str)
            
            if current_line:
                formatted_ranges.append(current_line)
            
            atom_list = " \\\n    ".join(formatted_ranges)
        
        return atom_list
    
    def _customize_plumed_outputs(self, plumed_content: str, replica_index: int) -> str:
        """
        Customize PLUMED output file names for replica
        
        Args:
            plumed_content: Original PLUMED content
            replica_index: Replica index
            
        Returns:
            Customized PLUMED content
        """
        # Replace common output file patterns with replica-specific names
        patterns = [
            (r'FILE=(\w+)\.(\w+)', rf'FILE=\1_replica{replica_index}.\2'),
            (r'STRIDE=(\d+)', r'STRIDE=\1'),  # Keep stride as is
        ]
        
        customized = plumed_content
        for pattern, replacement in patterns:
            customized = re.sub(pattern, replacement, customized)
        
        return customized
    
    def create_temperature_summary(self) -> None:
        """Create summary file with temperature information"""
        summary_path = Path(self.replica_data['base_output_dir']) / "temperature_summary.txt"
        
        summary_content = f"""# REST2 Temperature Summary
# Generated automatically

Total replicas: {self.n_replicas}
Temperature range: {min(self.temperatures):.1f} - {max(self.temperatures):.1f} K
Reference temperature: {self.T_ref:.1f} K
Scaling method: {self.replica_data.get('scaling_method', 'unknown')}

Replica Information:
{'Index':<6} {'Temperature (K)':<15} {'λ factor':<12} {'√λ factor':<12}
{'-'*50}
"""
        
        for replica in self.replicas:
            idx = replica['index']
            temp = replica['temperature']
            lambda_val = replica['scaling_factor']
            sqrt_lambda = np.sqrt(lambda_val)
            
            summary_content += f"{idx:<6} {temp:<15.1f} {lambda_val:<12.6f} {sqrt_lambda:<12.6f}\n"
        
        summary_content += f"""
{'-'*50}

REST2 Scaling Explanation:
- λ (lambda): Solute-solute interaction scaling factor
- √λ (sqrt lambda): Solute-solvent interaction scaling factor
- Solvent-solvent interactions remain unscaled (factor = 1.0)

Files generated for each replica:
- topol.top: Temperature-scaled topology file
- rest2.mdp: MDP file with replica-specific temperature
- plumed.dat: PLUMED file with replica-specific outputs (if applicable)
"""
        
        with open(summary_path, 'w') as f:
            f.write(summary_content)
    
    def validate_temperature_setup(self) -> bool:
        """
        Validate that temperature setup is complete (不再检查index.ndx)
        """
        errors = []
        for replica in self.replicas:
            replica_index = replica['index']
            input_dir = Path(replica['input_dir'])
            # 只检查topol.top、rest2.mdp、input.tpr
            required_files = {
                'topol.top': 'Temperature-scaled topology file',
                'rest2.mdp': 'MDP file with temperature settings',
                'input.tpr': 'GROMACS TPR file'
            }
            for filename, description in required_files.items():
                file_path = input_dir / filename
                if not file_path.exists():
                    errors.append(f"Replica {replica_index}: Missing {description} ({filename})")
            # 拓扑文件检查同前
            topology_file = input_dir / "topol.top"
            if topology_file.exists():
                pass
            else:
                errors.append(f"Replica {replica_index}: Missing topology file")
            # MDP温度检查同前
            mdp_file = input_dir / "rest2.mdp"
            if mdp_file.exists():
                with open(mdp_file, 'r') as f:
                    content = f.read()
                    expected_temp = f"{replica['temperature']:.1f}"
                    if f"ref_t                   = {expected_temp} {expected_temp}" not in content:
                        errors.append(f"Replica {replica_index}: MDP file missing correct temperature")
        if errors:
            print("Temperature setup validation errors:")
            for error in errors:
                print(f"  - {error}")
            return False
        print(f"Temperature setup validation passed for all {self.n_replicas} replicas")
        print("All replicas ready for REST2 simulation")
        return True
    
    def print_temperature_summary(self) -> None:
        """Print temperature setup summary"""
        print(f"Temperature setup:")
        print(f"  Replicas: {self.n_replicas}")
        print(f"  Temperature range: {min(self.temperatures):.1f} - {max(self.temperatures):.1f} K")
        print(f"  Reference temperature: {self.T_ref:.1f} K")


def main():
    """Test temperature controller"""
    try:
        # Mock config and replica data for testing
        class MockConfig:
            def __init__(self):
                self.params = {
                    'replex': 200,
                    'simulation.production_time': 100.0,
                    'md_results_dir': 'example/MD_results',
                    'plumed_dat': 'templates/plumed.dat'
                }
            
            def get_parameter(self, key, default=None):
                return self.params.get(key, default)
        
        mock_replica_data = {
            'replicas': [
                {'index': 0, 'temperature': 300.0, 'scaling_factor': 1.000000, 
                 'input_dir': './test_rest2/replica_0/input'},
                {'index': 1, 'temperature': 313.3, 'scaling_factor': 0.957447,
                 'input_dir': './test_rest2/replica_1/input'}
            ],
            'n_replicas': 2,
            'temperature_range': (300.0, 340.0),
            'scaling_method': 'linear',
            'base_output_dir': './test_rest2'
        }
        
        mock_solute_data = {
            'solute_atom_indices': [10, 11, 12, 15, 16, 17, 20, 21, 22, 25, 26, 27, 30],
            'target_atom_indices': [10, 11, 12, 15, 16, 17, 20, 21, 22, 25, 26, 27, 30],
            'nearby_residue_ids': [1, 2, 3]
        }
        
        # Initialize temperature controller
        config = MockConfig()
        controller = TemperatureController(config, mock_replica_data, mock_solute_data)
        
        # Print summary
        controller.print_temperature_summary()
        
        # Create temperature summary file
        controller.create_temperature_summary()
        
        print("Temperature controller module ready for use")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()