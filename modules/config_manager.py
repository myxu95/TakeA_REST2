#!/usr/bin/env python3
"""
REST2 Configuration Manager Module
Manages configuration parameters for REST2 enhanced sampling simulations
"""

import yaml
import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Union

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


class ConfigValidationError(Exception):
    """Configuration validation error"""
    pass


class ConfigManager:
    """
    REST2 Configuration Manager
    Handles loading, validation and access to configuration parameters
    """
    
    # Default configuration template
    DEFAULT_CONFIG = {
        # Basic REST2 settings
        'target_type': 'peptide',  # 'peptide' or 'small_molecule'
        'T_min': 300.0,            # Minimum temperature (K)
        'T_max': 340.0,            # Maximum temperature (K)
        'replex': 200,             # Exchange interval for gmx_mpi -replex
        'n_replicas': 8,           # Number of replicas
        'scaling_method': 'linear',  # 'exponential' or 'linear'
        
        # File paths for GROMACS input
        'input_tpr': None,         # -f: Input .tpr file (from completed MD)
        'topology': None,          # -p: Topology file (.top)
        'plumed_dat': 'templates/plumed.dat',  # PLUMED input file (.dat)
        'output_tpr': None,        # -o: Output topology file for REST2
        
        # Selection parameters
        'distance_range': 6.0,     # Cutoff distance for nearby residues (Angstrom)
        'target_selection': 'chain A',  # Target selection: 'chain A' for peptide or 'resname LIG' for ligand
        'use_trajectory': False,   # Whether to use trajectory for selection
        'occupancy_threshold': 0.5, # Threshold for trajectory-based selection
        
        # MD results directory structure
        'md_results_dir': 'example/MD_results',  # Base directory containing MD results
        
        # Optional parameters
        'output_dir': './rest2_simulation',
        'force_overwrite': False
    }
    
    # Required MD files structure template
    MD_FILES_STRUCTURE = {
        'required_files': [
            'md.tpr',           # Production MD tpr file
            'topol.top',        # System topology
            'md.gro'            # Final structure
        ],
        'optional_files': [
            'md.xtc',           # Trajectory (if use_trajectory=True)
            'md.edr',           # Energy file
            'md.log',           # Log file
            'index.ndx'         # Index file
        ]
    }
    
    def __init__(self, config_file: Optional[str] = None):
        """
        Initialize configuration manager
        
        Args:
            config_file: Path to configuration file
        """
        self.config_file = config_file
        self.config = {}
        
        if config_file:
            self.load_config(config_file)
        else:
            self.config = self.DEFAULT_CONFIG.copy()
    
    def load_config(self, config_file: str) -> None:
        """
        Load configuration from YAML file
        
        Args:
            config_file: Path to configuration file
        """
        config_path = Path(config_file)
        
        if not config_path.exists():
            self.create_default_config(config_file)
            raise FileNotFoundError(f"Config file created at {config_file}. Please edit and run again.")
        
        try:
            if FileUtils:
                # Use unified file utilities
                user_config = FileUtils.load_yaml(config_path)
            else:
                # Fallback to direct YAML loading
                with open(config_path, 'r', encoding='utf-8') as f:
                    user_config = yaml.safe_load(f) or {}
            
            # Merge with defaults
            self.config = self.DEFAULT_CONFIG.copy()
            self.config.update(user_config)
            
        except Exception as e:
            raise ConfigValidationError(f"Failed to load config: {e}")
    
    def create_default_config(self, config_file: str) -> None:
        """
        Create default configuration file
        
        Args:
            config_file: Path to configuration file
        """
        default_config = {
            'target_type': 'peptide',
            'T_min': 300.0,
            'T_max': 340.0,
            'replex': 200,
            'n_replicas': 8,
            'scaling_method': 'linear',
            'input_tpr': 'example/MD_results/md.tpr',
            'topology': 'example/MD_results/topol.top',
            'plumed_dat': 'templates/plumed.dat',
            'output_tpr': 'rest2_topol.top',
            'distance_range': 6.0,
            'target_selection': 'chain A',
            'use_trajectory': False,
            'occupancy_threshold': 0.5,
            'md_results_dir': 'example/MD_results',
            'output_dir': './rest2_simulation',
            'force_overwrite': False
        }
        
        if FileUtils:
            # Use unified file utilities
            FileUtils.save_yaml(default_config, config_file)
        else:
            # Fallback to direct YAML writing
            config_path = Path(config_file)
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with open(config_path, 'w', encoding='utf-8') as f:
                yaml.dump(default_config, f, default_flow_style=False, indent=2)
    
    def validate_config(self) -> None:
        """
        Validate configuration parameters using unified validation framework
        
        Raises:
            ConfigValidationError: If configuration is invalid
        """
        if ValidationFramework:
            # Use unified validation framework
            config_errors = ValidationFramework.validate_configuration(self.config)
            file_errors = ValidationFramework.validate_file_paths(self.config)
            
            all_errors = config_errors + file_errors
            
            if all_errors:
                error_msg = "\n".join(all_errors)
                raise ConfigValidationError(f"Configuration validation failed:\n{error_msg}")
        else:
            # Fallback validation
            errors = []
            
            # Basic parameter validation
            if self.config['T_min'] >= self.config['T_max']:
                errors.append("T_min must be less than T_max")
            
            if self.config['replex'] <= 0:
                errors.append("replex must be positive")
            
            if self.config['distance_range'] <= 0:
                errors.append("distance_range must be positive")
            
            if self.config['occupancy_threshold'] < 0 or self.config['occupancy_threshold'] > 1:
                errors.append("occupancy_threshold must be between 0 and 1")
            
            if errors:
                raise ConfigValidationError("\n".join(errors))
    
    def _auto_detect_files(self, md_dir: Path) -> None:
        """Auto-detect files in MD results directory"""
        if FileUtils:
            # Use unified file utilities
            file_mapping = {
                'input_tpr': ['md.tpr', 'topol.tpr', 'input.tpr'],
                'topology': ['topol.top', 'system.top', 'input.top'],
                'structure': ['md.gro', 'system.gro', 'input.gro']
            }
            
            detected_files = FileUtils.auto_detect_files(md_dir, file_mapping)
            
            for config_key, file_path in detected_files.items():
                if file_path and not self.config.get(config_key):
                    self.config[config_key] = file_path
        else:
            # Fallback auto-detection
            common_files = {
                'input_tpr': ['md.tpr', 'topol.tpr', 'input.tpr'],
                'topology': ['topol.top', 'system.top', 'input.top'],
                'structure': ['md.gro', 'system.gro', 'input.gro']
            }
            
            for config_key, possible_names in common_files.items():
                if not self.config.get(config_key):
                    for filename in possible_names:
                        file_path = md_dir / filename
                        if file_path.exists():
                            self.config[config_key] = str(file_path)
                            break
    
    def get_parameter(self, key: str, default: Any = None) -> Any:
        """
        Get configuration parameter
        
        Args:
            key: Parameter key
            default: Default value if key not found
            
        Returns:
            Parameter value
        """
        return self.config.get(key, default)
    
    def set_parameter(self, key: str, value: Any) -> None:
        """
        Set configuration parameter
        
        Args:
            key: Parameter key
            value: Parameter value
        """
        self.config[key] = value
    
    def get_all_parameters(self) -> Dict[str, Any]:
        """Get all configuration parameters"""
        return self.config.copy()
    
    def get_temperature_ladder(self) -> list:
        """
        Get temperature ladder using unified calculator
        
        Returns:
            List of temperatures for each replica
        """
        if TemperatureCalculator:
            return TemperatureCalculator.calculate_temperature_ladder(
                self.config['T_min'],
                self.config['T_max'],
                self.config['n_replicas'],
                self.config['scaling_method']
            )
        else:
            # Fallback calculation
            T_min = self.config['T_min']
            T_max = self.config['T_max']
            n_replicas = self.config['n_replicas']
            
            if n_replicas == 1:
                return [T_min]
            
            step = (T_max - T_min) / (n_replicas - 1)
            return [T_min + i * step for i in range(n_replicas)]
    
    def get_scaling_factors(self) -> list:
        """
        Get REST2 scaling factors using unified calculator
        
        Returns:
            List of scaling factors for each replica
        """
        temperatures = self.get_temperature_ladder()
        
        if TemperatureCalculator:
            return TemperatureCalculator.calculate_scaling_factors(temperatures)
        else:
            # Fallback calculation
            T_ref = temperatures[0]
            return [T_ref / T for T in temperatures]
    
    def print_summary(self) -> None:
        """Print configuration summary using unified formatter"""
        if OutputFormatter:
            # Use unified output formatter
            OutputFormatter.print_configuration_summary(self.config)
            
            # Print temperature ladder using unified calculator
            temperatures = self.get_temperature_ladder()
            scaling_factors = self.get_scaling_factors()
            
            if TemperatureCalculator:
                TemperatureCalculator.print_temperature_summary(
                    temperatures, scaling_factors, self.config['scaling_method']
                )
            else:
                OutputFormatter.print_temperature_summary(
                    temperatures, scaling_factors, self.config['scaling_method']
                )
        else:
            # Fallback summary
            print("\n" + "="*50)
            print("REST2 Configuration Summary")
            print("="*50)
            
            print(f"Target type: {self.config['target_type']}")
            print(f"Target selection: {self.config['target_selection']}")
            print(f"Temperature range: {self.config['T_min']:.1f} - {self.config['T_max']:.1f} K")
            print(f"Number of replicas: {self.config['n_replicas']}")
            print(f"Exchange interval: {self.config['replex']} steps")
            print(f"Scaling method: {self.config['scaling_method']}")
            print(f"Distance cutoff: {self.config['distance_range']} Angstrom")
            
            if self.config['use_trajectory']:
                print(f"Using trajectory-based selection (threshold: {self.config['occupancy_threshold']})")
            else:
                print("Using static structure for selection")
            
            # Print temperature ladder
            temperatures = self.get_temperature_ladder()
            scaling_factors = self.get_scaling_factors()
            
            print("\nReplica temperatures and scaling factors:")
            for i, (T, lambda_val) in enumerate(zip(temperatures, scaling_factors)):
                print(f"  Replica {i}: T={T:.1f}K, λ={lambda_val:.3f}")
            
            print(f"\nInput files:")
            print(f"  TPR file: {self.config.get('input_tpr', 'Not specified')}")
            print(f"  Topology: {self.config.get('topology', 'Not specified')}")
            print(f"  PLUMED: {self.config.get('plumed_dat', 'Not specified')}")
            
            print(f"Output directory: {self.config['output_dir']}")
            print("="*50 + "\n")
    
    def get_md_files_info(self) -> Dict[str, Any]:
        """
        Get information about expected MD files structure
        
        Returns:
            Dictionary with file structure information
        """
        return self.MD_FILES_STRUCTURE.copy()


def main():
    """Test configuration manager"""
    try:
        # Test with default config
        config_manager = ConfigManager()
        
        print("Testing Configuration Manager")
        print("=" * 40)
        
        # Test temperature calculations
        print("Temperature calculations:")
        temperatures = config_manager.get_temperature_ladder()
        scaling_factors = config_manager.get_scaling_factors()
        
        print(f"Temperatures: {temperatures}")
        print(f"Scaling factors: {scaling_factors}")
        
        # Test parameter access
        print(f"\nT_min: {config_manager.get_parameter('T_min')}")
        print(f"n_replicas: {config_manager.get_parameter('n_replicas')}")
        
        # Test validation
        print("\nTesting validation:")
        try:
            config_manager.validate_config()
            print("✓ Configuration validation passed")
        except ConfigValidationError as e:
            print(f"✗ Configuration validation failed: {e}")
        
        # Test summary
        config_manager.print_summary()
        
        print("✓ Configuration manager tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()