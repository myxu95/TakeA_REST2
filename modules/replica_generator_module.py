#!/usr/bin/env python3
"""
REST2 Replica Generator Module
Generates replica directories and basic files for REST2 enhanced sampling
"""

import os
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional
import numpy as np


class ReplicaGeneratorError(Exception):
    """Replica generator error"""
    pass


class ReplicaGenerator:
    """
    Replica generator for REST2 simulations
    Creates temperature ladder and replica directory structure
    """
    
    def __init__(self, config_manager):
        """
        Initialize replica generator
        
        Args:
            config_manager: ConfigManager instance with simulation parameters
        """
        self.config = config_manager
        self.T_min = config_manager.get_parameter('T_min')
        self.T_max = config_manager.get_parameter('T_max')
        self.n_replicas = config_manager.get_parameter('n_replicas')
        self.scaling_method = config_manager.get_parameter('scaling_method')
        self.output_dir = Path(config_manager.get_parameter('output_dir'))
        
        # Temperature ladder and scaling factors
        self.temperatures = []
        self.scaling_factors = []
        
        # Generate temperature ladder
        self._calculate_temperature_ladder()
    
    def _calculate_temperature_ladder(self) -> None:
        """Calculate temperature ladder for replicas"""
        if self.n_replicas == 1:
            self.temperatures = [self.T_min]
        elif self.scaling_method == 'linear':
            # Linear temperature spacing
            self.temperatures = np.linspace(self.T_min, self.T_max, self.n_replicas).tolist()
        else:  # exponential
            # Exponential temperature spacing
            if self.T_min <= 0:
                raise ReplicaGeneratorError("T_min must be positive for exponential scaling")
            
            ratio = (self.T_max / self.T_min) ** (1.0 / (self.n_replicas - 1))
            self.temperatures = [self.T_min * (ratio ** i) for i in range(self.n_replicas)]
        
        # Calculate REST2 scaling factors (λ = T_ref/T)
        T_ref = self.temperatures[0]  # Reference temperature (lowest)
        self.scaling_factors = [T_ref / T for T in self.temperatures]
    
    def get_temperature_ladder(self) -> List[float]:
        """Get temperature ladder"""
        return self.temperatures.copy()
    
    def get_scaling_factors(self) -> List[float]:
        """Get REST2 scaling factors"""
        return self.scaling_factors.copy()
    
    def setup_replica_directories(self) -> None:
        """Create directory structure for all replicas"""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        for i in range(self.n_replicas):
            replica_dir = self.output_dir / f"replica_{i}"
            replica_dir.mkdir(exist_ok=True)
            
            # Create subdirectories for organization
            (replica_dir / "input").mkdir(exist_ok=True)
            (replica_dir / "output").mkdir(exist_ok=True)
    
    def copy_base_files_to_replicas(self, modified_topology: str) -> None:
        """
        Copy base input files to each replica directory
        
        Args:
            modified_topology: Path to REST2-modified topology file
        """
        input_tpr = self.config.get_parameter('input_tpr')
        plumed_dat = self.config.get_parameter('plumed_dat')
        
        # Validate input files
        required_files = {
            'input_tpr': input_tpr,
            'modified_topology': modified_topology
        }
        
        for file_type, file_path in required_files.items():
            if not file_path or not Path(file_path).exists():
                raise FileNotFoundError(f"Required file not found: {file_path} ({file_type})")
        
        # Copy files to each replica
        for i in range(self.n_replicas):
            replica_input_dir = self.output_dir / f"replica_{i}" / "input"
            
            # Copy TPR file
            shutil.copy2(input_tpr, replica_input_dir / "input.tpr")
            
            # Copy modified topology
            shutil.copy2(modified_topology, replica_input_dir / "topol.top")
            
            # Copy PLUMED file if provided
            if plumed_dat and Path(plumed_dat).exists():
                shutil.copy2(plumed_dat, replica_input_dir / "plumed.dat")
    
    def create_replica_info_files(self) -> None:
        """Create information files for each replica"""
        for i in range(self.n_replicas):
            replica_dir = self.output_dir / f"replica_{i}"
            
            # Create replica info file
            info_content = f"""# Replica {i} Information
# Generated automatically for REST2 simulation

Replica Index: {i}
Temperature: {self.temperatures[i]:.1f} K
REST2 Scaling Factor (λ): {self.scaling_factors[i]:.6f}

# Directory Structure:
#   input/     - Input files (TPR, topology, PLUMED)
#   output/    - Output files (will be created during simulation)

# Next Steps:
# 1. Temperature controller will modify input files
# 2. MDP files will be generated with appropriate settings
# 3. Run scripts will be created for job submission
"""
            
            info_file = replica_dir / "replica_info.txt"
            with open(info_file, 'w') as f:
                f.write(info_content)
    
    def get_replica_data(self) -> Dict[str, Any]:
        """
        Get replica data for use by other modules
        
        Returns:
            Dictionary with replica information
        """
        replica_data = []
        
        for i in range(self.n_replicas):
            replica_info = {
                'index': i,
                'temperature': self.temperatures[i],
                'scaling_factor': self.scaling_factors[i],
                'directory': str(self.output_dir / f"replica_{i}"),
                'input_dir': str(self.output_dir / f"replica_{i}" / "input"),
                'output_dir': str(self.output_dir / f"replica_{i}" / "output")
            }
            replica_data.append(replica_info)
        
        return {
            'replicas': replica_data,
            'n_replicas': self.n_replicas,
            'temperature_range': (self.T_min, self.T_max),
            'scaling_method': self.scaling_method,
            'base_output_dir': str(self.output_dir)
        }
    
    def print_replica_summary(self) -> None:
        """Print summary of replica setup"""
        print("\n" + "="*70)
        print("REST2 Replica Configuration Summary")
        print("="*70)
        
        print(f"Number of replicas: {self.n_replicas}")
        print(f"Temperature range: {self.T_min:.1f} - {self.T_max:.1f} K")
        print(f"Scaling method: {self.scaling_method}")
        print()
        
        print("Replica Temperature Ladder:")
        print(f"{'Replica':<8} {'Temperature (K)':<15} {'λ factor':<12}")
        print("-"*40)
        
        for i in range(self.n_replicas):
            print(f"{i:<8} {self.temperatures[i]:<15.1f} {self.scaling_factors[i]:<12.6f}")
        
        print()
        print(f"Base output directory: {self.output_dir}")
        print("="*70)
        print("Directory structure created for each replica:")
        print("  replica_X/")
        print("    ├── input/          # Input files (TPR, topology, PLUMED)")
        print("    ├── output/         # Output files (created during simulation)")
        print("    └── replica_info.txt # Replica information")
        print("="*70 + "\n")
    
    def validate_replica_setup(self) -> bool:
        """
        Validate that replica setup is complete
        
        Returns:
            True if validation passes
        """
        errors = []
        
        # Check main output directory
        if not self.output_dir.exists():
            errors.append(f"Output directory does not exist: {self.output_dir}")
        
        # Check each replica directory
        for i in range(self.n_replicas):
            replica_dir = self.output_dir / f"replica_{i}"
            input_dir = replica_dir / "input"
            output_dir = replica_dir / "output"
            info_file = replica_dir / "replica_info.txt"
            
            if not replica_dir.exists():
                errors.append(f"Replica directory missing: replica_{i}")
                continue
            
            if not input_dir.exists():
                errors.append(f"Input directory missing: replica_{i}/input")
            
            if not output_dir.exists():
                errors.append(f"Output directory missing: replica_{i}/output")
            
            if not info_file.exists():
                errors.append(f"Info file missing: replica_{i}/replica_info.txt")
            
            # Check basic input files (TPR and topology should be copied)
            basic_files = ["input.tpr", "topol.top"]
            for filename in basic_files:
                file_path = input_dir / filename
                if not file_path.exists():
                    errors.append(f"Missing basic file: replica_{i}/input/{filename}")
        
        if errors:
            print("Validation errors:")
            for error in errors:
                print(f"  - {error}")
            return False
        
        print(f"Validation passed: All {self.n_replicas} replicas properly set up")
        print("Ready for temperature controller and MDP generation")
        return True
    
    def get_multidir_string(self) -> str:
        """
        Get multidir string for GROMACS commands
        
        Returns:
            Space-separated string of replica directories
        """
        return " ".join([f"replica_{i}" for i in range(self.n_replicas)])


def main():
    """Test replica generator"""
    try:
        # Mock config manager for testing
        class MockConfig:
            def __init__(self):
                self.params = {
                    'T_min': 300.0,
                    'T_max': 340.0,
                    'n_replicas': 4,
                    'scaling_method': 'linear',
                    'replex': 200,
                    'output_dir': './test_rest2',
                    'input_tpr': 'example/MD_results/md.tpr',
                    'plumed_dat': 'templates/plumed.dat'
                }
            
            def get_parameter(self, key, default=None):
                return self.params.get(key, default)
        
        # Initialize replica generator
        config = MockConfig()
        generator = ReplicaGenerator(config)
        
        # Print summary
        generator.print_replica_summary()
        
        # Set up directories
        generator.setup_replica_directories()
        
        # Create info files
        generator.create_replica_info_files()
        
        # Get replica data for other modules
        replica_data = generator.get_replica_data()
        print(f"Replica data prepared for next modules:")
        print(f"  - {replica_data['n_replicas']} replicas configured")
        print(f"  - Temperature range: {replica_data['temperature_range']}")
        print(f"  - Multidir string: {generator.get_multidir_string()}")
        
        # Validate setup
        generator.validate_replica_setup()
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()