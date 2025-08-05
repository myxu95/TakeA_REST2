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
import sys

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
        """Calculate temperature ladder using unified calculator"""
        try:
            # Use unified temperature calculator
            self.temperatures, self.scaling_factors = TemperatureCalculator.calculate_temperature_and_scaling(
                self.T_min, self.T_max, self.n_replicas, self.scaling_method
            )
        except TemperatureCalculationError as e:
            raise ReplicaGeneratorError(f"Temperature calculation failed: {e}")
    
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
            
            shutil.copy2(input_tpr, replica_input_dir / "input.tpr")
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
Replica Index: {i}
Temperature: {self.temperatures[i]:.1f} K
REST2 Scaling Factor (λ): {self.scaling_factors[i]:.6f}
"""
            
            info_file = replica_dir / "replica_info.txt"
            with open(info_file, 'w') as f:
                f.write(info_content)
    
    def get_replica_data(self) -> Dict[str, Any]:
        """
        Get replica data for other modules
        
        Returns:
            Dictionary containing replica information
        """
        replicas = []
        for i in range(self.n_replicas):
            replica_data = {
                'index': i,
                'temperature': self.temperatures[i],
                'scaling_factor': self.scaling_factors[i],
                'input_dir': str(self.output_dir / f"replica_{i}" / "input"),
                'output_dir': str(self.output_dir / f"replica_{i}" / "output"),
                'replica_dir': str(self.output_dir / f"replica_{i}")
            }
            replicas.append(replica_data)
        
        return {
            'n_replicas': self.n_replicas,
            'temperatures': self.temperatures,
            'scaling_factors': self.scaling_factors,
            'replicas': replicas,
            'base_output_dir': str(self.output_dir)
        }
    
    def print_replica_summary(self) -> None:
        """Print replica generation summary using unified formatter"""
        if OutputFormatter:
            # Use unified output formatter
            replica_data = self.get_replica_data()
            OutputFormatter.print_replica_summary(replica_data)
            
            # Print temperature summary
            if TemperatureCalculator:
                TemperatureCalculator.print_temperature_summary(
                    self.temperatures, self.scaling_factors, self.scaling_method
                )
            else:
                OutputFormatter.print_temperature_summary(
                    self.temperatures, self.scaling_factors, self.scaling_method
                )
        else:
            # Fallback summary
            print("\n" + "="*50)
            print("Replica Generation Summary")
            print("="*50)
            
            print(f"Number of replicas: {self.n_replicas}")
            print(f"Output directory: {self.output_dir}")
            print(f"Scaling method: {self.scaling_method}")
            
            # Print temperature ladder
            print("\nReplica temperatures and scaling factors:")
            for i, (T, lambda_val) in enumerate(zip(self.temperatures, self.scaling_factors)):
                print(f"  Replica {i}: T={T:.1f}K, λ={lambda_val:.6f}")
            
            print("\nDirectory structure:")
            for i in range(self.n_replicas):
                replica_dir = self.output_dir / f"replica_{i}"
                print(f"  replica_{i}/")
                print(f"    ├── input/     (TPR, topology, PLUMED)")
                print(f"    ├── output/    (simulation results)")
                print(f"    └── replica_info.txt")
            
            print("="*50)
    
    def validate_replica_setup(self) -> bool:
        """
        Validate replica setup using unified validation framework
        
        Returns:
            True if setup is valid
        """
        try:
            if ValidationFramework:
                # Use unified validation framework
                errors = ValidationFramework.validate_replica_setup(
                    self.output_dir, self.n_replicas, self.temperatures, self.scaling_factors
                )
                
                ValidationFramework.print_validation_summary(errors, "Replica setup validation")
                return len(errors) == 0
            else:
                # Fallback validation
                # Check output directory exists
                if not self.output_dir.exists():
                    print(f"✗ Output directory does not exist: {self.output_dir}")
                    return False
                
                # Check all replica directories exist
                for i in range(self.n_replicas):
                    replica_dir = self.output_dir / f"replica_{i}"
                    if not replica_dir.exists():
                        print(f"✗ Replica directory does not exist: {replica_dir}")
                        return False
                
                print("✓ Replica setup validation passed")
                return True
                
        except Exception as e:
            print(f"✗ Replica setup validation failed: {e}")
            return False
    
    def get_multidir_string(self) -> str:
        """
        Get multidir string for GROMACS
        
        Returns:
            Space-separated list of replica directories
        """
        return " ".join([f"replica_{i}" for i in range(self.n_replicas)])


def main():
    """Test replica generator"""
    try:
        # Create mock config for testing
        class MockConfig:
            def __init__(self):
                self.params = {
                    'T_min': 300.0,
                    'T_max': 340.0,
                    'n_replicas': 8,
                    'scaling_method': 'linear',
                    'output_dir': './test_replicas',
                    'input_tpr': 'test_input.tpr',
                    'plumed_dat': 'test_plumed.dat'
                }
            
            def get_parameter(self, key, default=None):
                return self.params.get(key, default)
        
        # Test replica generator
        config = MockConfig()
        generator = ReplicaGenerator(config)
        
        print("Testing Replica Generator")
        print("=" * 40)
        
        # Test temperature calculations
        temperatures = generator.get_temperature_ladder()
        scaling_factors = generator.get_scaling_factors()
        
        print(f"Temperatures: {temperatures}")
        print(f"Scaling factors: {scaling_factors}")
        
        # Test summary
        generator.print_replica_summary()
        
        # Test validation (will fail since directories don't exist)
        print("\nTesting validation (expected to fail):")
        generator.validate_replica_setup()
        
        print("✓ Replica generator tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()