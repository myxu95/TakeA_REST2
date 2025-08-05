#!/usr/bin/env python3
"""
REST2 Validation Framework Module
Unified validation framework for REST2 enhanced sampling simulations
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Union

# Import from same package
from .temperature_calculator import TemperatureCalculator, TemperatureCalculationError


class ValidationError(Exception):
    """Validation error"""
    pass


class ValidationFramework:
    """
    Unified validation framework for REST2 simulations
    Handles all validation logic across different modules
    """
    
    @staticmethod
    def validate_configuration(config: Dict[str, Any]) -> List[str]:
        """
        Validate configuration parameters
        
        Args:
            config: Configuration dictionary
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Validate temperature parameters
        if TemperatureCalculator:
            try:
                TemperatureCalculator.validate_temperature_parameters(
                    config.get('T_min', 0),
                    config.get('T_max', 0),
                    config.get('n_replicas', 0),
                    config.get('scaling_method', 'linear')
                )
            except TemperatureCalculationError as e:
                errors.append(f"Temperature parameters: {e}")
        
        # Validate other parameters
        if config.get('replex', 0) <= 0:
            errors.append("replex must be positive")
        
        if config.get('distance_range', 0) <= 0:
            errors.append("distance_range must be positive")
        
        occupancy_threshold = config.get('occupancy_threshold', 0.5)
        if occupancy_threshold < 0 or occupancy_threshold > 1:
            errors.append("occupancy_threshold must be between 0 and 1")
        
        # Validate target selection
        target_selection = config.get('target_selection', '')
        if not target_selection:
            errors.append("target_selection must be specified")
        
        # Validate target type
        target_type = config.get('target_type', '')
        if target_type not in ['peptide', 'small_molecule']:
            errors.append("target_type must be 'peptide' or 'small_molecule'")
        
        return errors
    
    @staticmethod
    def validate_file_paths(config: Dict[str, Any]) -> List[str]:
        """
        Validate file paths in configuration
        
        Args:
            config: Configuration dictionary
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check required files if they are specified
        required_files = ['input_tpr', 'topology']
        
        for file_key in required_files:
            file_path = config.get(file_key)
            if file_path and not Path(file_path).exists():
                errors.append(f"File not found: {file_path} ({file_key})")
        
        # Check MD results directory
        md_dir = Path(config.get('md_results_dir', ''))
        if not md_dir.exists():
            errors.append(f"MD results directory not found: {md_dir}")
        
        # Check trajectory file if needed
        if config.get('use_trajectory', False):
            traj_file = md_dir / 'md.xtc'
            if not traj_file.exists():
                errors.append("Trajectory file (md.xtc) required but not found")
        
        return errors
    
    @staticmethod
    def validate_structure_files(structure_file: str, topology_file: str, 
                               trajectory_file: Optional[str] = None) -> List[str]:
        """
        Validate structure analysis input files
        
        Args:
            structure_file: Structure file path
            topology_file: Topology file path
            trajectory_file: Trajectory file path (optional)
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check structure file
        if not Path(structure_file).exists():
            errors.append(f"Structure file not found: {structure_file}")
        
        # Check topology file
        if not Path(topology_file).exists():
            errors.append(f"Topology file not found: {topology_file}")
        
        # Check trajectory file if provided
        if trajectory_file and not Path(trajectory_file).exists():
            errors.append(f"Trajectory file not found: {trajectory_file}")
        
        return errors
    
    @staticmethod
    def validate_topology_modification(input_topology: str, output_topology: str,
                                    solute_data: Dict[str, Any]) -> List[str]:
        """
        Validate topology modification
        
        Args:
            input_topology: Input topology file path
            output_topology: Output topology file path
            solute_data: Solute selection data
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check input topology exists
        if input_topology and not Path(input_topology).exists():
            errors.append(f"Input topology file not found: {input_topology}")
        
        # Check output topology was created
        if not Path(output_topology).exists():
            errors.append(f"Output topology file not created: {output_topology}")
        
        # Validate solute data
        required_keys = ['target_atom_indices', 'nearby_residue_ids', 'solute_atom_indices']
        for key in required_keys:
            if key not in solute_data:
                errors.append(f"Missing solute data key: {key}")
        
        # Check solute data is not empty
        if solute_data.get('solute_atom_indices'):
            if len(solute_data['solute_atom_indices']) == 0:
                errors.append("No solute atoms selected")
        
        return errors
    
    @staticmethod
    def validate_replica_setup(output_dir: Path, n_replicas: int, 
                             temperatures: List[float], scaling_factors: List[float]) -> List[str]:
        """
        Validate replica setup
        
        Args:
            output_dir: Output directory path
            n_replicas: Number of replicas
            temperatures: List of temperatures
            scaling_factors: List of scaling factors
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check output directory exists
        if not output_dir.exists():
            errors.append(f"Output directory does not exist: {output_dir}")
            return errors
        
        # Check all replica directories exist
        for i in range(n_replicas):
            replica_dir = output_dir / f"replica_{i}"
            if not replica_dir.exists():
                errors.append(f"Replica directory does not exist: {replica_dir}")
                continue
            
            # Check subdirectories
            input_dir = replica_dir / "input"
            output_dir_replica = replica_dir / "output"
            
            if not input_dir.exists():
                errors.append(f"Input directory does not exist: {input_dir}")
            
            if not output_dir_replica.exists():
                errors.append(f"Output directory does not exist: {output_dir_replica}")
            
            # Check required files in input directory
            required_files = ['input.tpr', 'topol.top']
            for filename in required_files:
                file_path = input_dir / filename
                if not file_path.exists():
                    errors.append(f"Required file not found: {file_path}")
        
        # Check temperature calculations
        if len(temperatures) != n_replicas:
            errors.append(f"Temperature count mismatch: {len(temperatures)} != {n_replicas}")
        
        if len(scaling_factors) != n_replicas:
            errors.append(f"Scaling factor count mismatch: {len(scaling_factors)} != {n_replicas}")
        
        return errors
    
    @staticmethod
    def validate_temperature_setup(replica_data: Dict[str, Any], 
                                 solute_data: Optional[Dict[str, Any]] = None) -> List[str]:
        """
        Validate temperature control setup
        
        Args:
            replica_data: Replica data dictionary
            solute_data: Solute data dictionary (optional)
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        replicas = replica_data.get('replicas', [])
        n_replicas = replica_data.get('n_replicas', 0)
        
        if len(replicas) != n_replicas:
            errors.append(f"Replica count mismatch: {len(replicas)} != {n_replicas}")
        
        # Check each replica has required files
        for i, replica in enumerate(replicas):
            input_dir = Path(replica.get('input_dir', ''))
            
            if not input_dir.exists():
                errors.append(f"Replica {i} input directory not found: {input_dir}")
                continue
            
            # Check required files
            required_files = ['input.tpr', 'topol.top']
            for filename in required_files:
                file_path = input_dir / filename
                if not file_path.exists():
                    errors.append(f"Replica {i} missing file: {file_path}")
        
        # Validate solute data if provided
        if solute_data:
            if 'solute_atom_indices' not in solute_data:
                errors.append("Missing solute atom indices in solute data")
            elif len(solute_data['solute_atom_indices']) == 0:
                errors.append("No solute atoms defined")
        
        return errors
    
    @staticmethod
    def validate_script_generation(config: Dict[str, Any], replica_data: Dict[str, Any]) -> List[str]:
        """
        Validate script generation setup
        
        Args:
            config: Configuration dictionary
            replica_data: Replica data dictionary
            
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check GROMACS settings
        gmx_command = config.get('gromacs', {}).get('gmx_mpi_command', 'gmx_mpi')
        if not gmx_command:
            errors.append("GROMACS command not specified")
        
        # Check replica data
        n_replicas = replica_data.get('n_replicas', 0)
        if n_replicas <= 0:
            errors.append("Invalid number of replicas")
        
        # Check output directory
        base_output_dir = Path(replica_data.get('base_output_dir', ''))
        if not base_output_dir.exists():
            errors.append(f"Base output directory not found: {base_output_dir}")
        
        return errors
    
    @staticmethod
    def print_validation_summary(errors: List[str], context: str = "Validation") -> None:
        """
        Print validation summary
        
        Args:
            errors: List of validation errors
            context: Validation context
        """
        if not errors:
            print(f"✓ {context} passed")
        else:
            print(f"✗ {context} failed:")
            for error in errors:
                print(f"  - {error}")
    
    @staticmethod
    def validate_complete_setup(config: Dict[str, Any], replica_data: Optional[Dict[str, Any]] = None,
                              solute_data: Optional[Dict[str, Any]] = None) -> bool:
        """
        Perform complete validation of REST2 setup
        
        Args:
            config: Configuration dictionary
            replica_data: Replica data dictionary (optional)
            solute_data: Solute data dictionary (optional)
            
        Returns:
            True if all validations pass
        """
        all_errors = []
        
        # Validate configuration
        config_errors = ValidationFramework.validate_configuration(config)
        all_errors.extend(config_errors)
        ValidationFramework.print_validation_summary(config_errors, "Configuration validation")
        
        # Validate file paths
        file_errors = ValidationFramework.validate_file_paths(config)
        all_errors.extend(file_errors)
        ValidationFramework.print_validation_summary(file_errors, "File path validation")
        
        # Validate replica setup if provided
        if replica_data:
            output_dir = Path(replica_data.get('base_output_dir', ''))
            n_replicas = replica_data.get('n_replicas', 0)
            temperatures = replica_data.get('temperatures', [])
            scaling_factors = replica_data.get('scaling_factors', [])
            
            replica_errors = ValidationFramework.validate_replica_setup(
                output_dir, n_replicas, temperatures, scaling_factors
            )
            all_errors.extend(replica_errors)
            ValidationFramework.print_validation_summary(replica_errors, "Replica setup validation")
        
        # Validate temperature setup if provided
        if replica_data and solute_data:
            temp_errors = ValidationFramework.validate_temperature_setup(replica_data, solute_data)
            all_errors.extend(temp_errors)
            ValidationFramework.print_validation_summary(temp_errors, "Temperature setup validation")
        
        # Print overall result
        if all_errors:
            print(f"\n✗ Complete validation failed with {len(all_errors)} errors")
            return False
        else:
            print(f"\n✓ Complete validation passed")
            return True


def main():
    """Test validation framework"""
    try:
        print("Testing Validation Framework")
        print("=" * 40)
        
        # Test configuration validation
        test_config = {
            'T_min': 300.0,
            'T_max': 340.0,
            'n_replicas': 8,
            'scaling_method': 'linear',
            'replex': 200,
            'distance_range': 6.0,
            'occupancy_threshold': 0.5,
            'target_selection': 'chain A',
            'target_type': 'peptide',
            'md_results_dir': 'example/MD_results'
        }
        
        print("Testing configuration validation:")
        errors = ValidationFramework.validate_configuration(test_config)
        ValidationFramework.print_validation_summary(errors, "Configuration validation")
        
        # Test temperature validation
        if TemperatureCalculator:
            print("\nTesting temperature validation:")
            try:
                TemperatureCalculator.validate_temperature_parameters(300.0, 340.0, 8, 'linear')
                print("✓ Temperature validation passed")
            except TemperatureCalculationError as e:
                print(f"✗ Temperature validation failed: {e}")
        
        # Test complete validation
        print("\nTesting complete validation:")
        ValidationFramework.validate_complete_setup(test_config)
        
        print("\n✓ Validation framework tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 