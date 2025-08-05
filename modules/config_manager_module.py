#!/usr/bin/env python3
"""
REST2 Configuration Manager Module
Manages configuration parameters for REST2 enhanced sampling simulations
"""

import yaml
import os
from pathlib import Path
from typing import Dict, Any, Optional, Union


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
            with open(config_path, 'r', encoding='utf-8') as f:
                user_config = yaml.safe_load(f) or {}
            
            # Merge with defaults
            self.config = {**self.DEFAULT_CONFIG, **user_config}
            self.config_file = str(config_path)
            
            # Validate configuration
            self.validate_config()
            
        except yaml.YAMLError as e:
            raise ConfigValidationError(f"YAML parsing error: {e}")
        except Exception as e:
            raise ConfigValidationError(f"Config loading failed: {e}")
    
    def create_default_config(self, config_file: str) -> None:
        """
        Create default configuration file with comments
        
        Args:
            config_file: Path for new configuration file
        """
        config_path = Path(config_file)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        config_content = """# REST2 Enhanced Sampling Configuration

# Target molecule type: 'peptide' or 'small_molecule'
target_type: peptide

# Temperature settings (K)
T_min: 300.0
T_max: 340.0

# Replica exchange settings
replex: 200         # Exchange interval (steps)
n_replicas: 8       # Number of replicas
scaling_method: linear  # 'exponential' or 'linear'

# GROMACS input files (will be auto-set from md_results_dir)
input_tpr: null     # Will be set to md_results_dir/md.tpr
topology: null      # Will be set to md_results_dir/topol.top
plumed_dat: templates/plumed.dat
output_tpr: ./rest2_system.tpr

# GROMACS execution settings
gromacs:
  gmx_mpi_command: gmx_mpi
  n_cpus: null      # Default: same as n_replicas
  n_gpus: null      # Default: same as n_replicas
  script_types:     # Scripts to generate
    - slurm
    - localrun
    - test

# Selection parameters
distance_range: 6.0         # Cutoff distance (Angstrom)
target_selection: chain A   # Target selection: 'chain A' for peptide or 'resname LIG' for ligand
use_trajectory: false       # Use trajectory for dynamic selection
occupancy_threshold: 0.5    # Occupancy threshold (if use_trajectory=true)

# MD results directory (containing completed MD simulation)
md_results_dir: example/MD_results

# Output settings
output_dir: ./rest2_simulation
force_overwrite: false

# Expected MD directory structure:
# md_results_dir/
# ├── md.tpr          # Production MD tpr file (required)
# ├── topol.top       # System topology (required)
# ├── md.gro          # Final structure (required)
# ├── md.xtc          # Trajectory (optional, needed if use_trajectory=true)
# ├── md.edr          # Energy file (optional)
# ├── md.log          # Log file (optional)
# └── index.ndx       # Index file (optional)
"""
        
        with open(config_path, 'w', encoding='utf-8') as f:
            f.write(config_content)
    
    def validate_config(self) -> None:
        """Validate configuration parameters"""
        errors = []
        
        # Basic parameter validation
        if self.config['target_type'] not in ['peptide', 'small_molecule']:
            errors.append("target_type must be 'peptide' or 'small_molecule'")
        
        if self.config['T_min'] >= self.config['T_max']:
            errors.append("T_min must be less than T_max")
        
        if self.config['T_min'] <= 0 or self.config['T_max'] <= 0:
            errors.append("Temperatures must be positive")
        
        if self.config['n_replicas'] < 2:
            errors.append("n_replicas must be at least 2")
        
        if self.config['replex'] < 1:
            errors.append("replex must be positive")
        
        if self.config['scaling_method'] not in ['exponential', 'linear']:
            errors.append("scaling_method must be 'exponential' or 'linear'")
        
        if self.config['distance_range'] <= 0:
            errors.append("distance_range must be positive")
        
        if not self.config['target_selection']:
            errors.append("target_selection must be specified")
        
        if self.config['use_trajectory'] and not (0 < self.config['occupancy_threshold'] <= 1):
            errors.append("occupancy_threshold must be between 0 and 1")
        
        # File existence validation
        self._validate_file_paths(errors)
        
        if errors:
            error_msg = "\n".join(f"  - {error}" for error in errors)
            raise ConfigValidationError(f"Configuration validation failed:\n{error_msg}")
    
    def _validate_file_paths(self, errors: list) -> None:
        """Validate required file paths"""
        # Check MD results directory is specified
        if not self.config['md_results_dir']:
            errors.append("md_results_dir must be specified")
            return
            
        md_dir = Path(self.config['md_results_dir'])
        if not md_dir.exists():
            errors.append(f"MD results directory does not exist: {md_dir}")
            return
        
        # Set file paths based on MD results directory
        self.config['input_tpr'] = str(md_dir / 'md.tpr')
        self.config['topology'] = str(md_dir / 'topol.top')
        
        # Validate required files exist
        required_files = ['input_tpr', 'topology']
        for file_key in required_files:
            file_path = self.config[file_key]
            if not Path(file_path).exists():
                errors.append(f"Required file does not exist: {file_path}")
        
        # Validate trajectory file if needed
        if self.config['use_trajectory']:
            traj_file = md_dir / 'md.xtc'
            if not traj_file.exists():
                errors.append("Trajectory file (md.xtc) required but not found")
    
    def _auto_detect_files(self, md_dir: Path) -> None:
        """Auto-detect file paths from MD results directory"""
        # This function is no longer used
        pass
    
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
        Calculate temperature ladder for replicas
        
        Returns:
            List of temperatures for each replica
        """
        T_min = self.config['T_min']
        T_max = self.config['T_max']
        n_replicas = self.config['n_replicas']
        method = self.config['scaling_method']
        
        if n_replicas == 1:
            return [T_min]
        
        if method == 'linear':
            # Linear spacing
            temperatures = []
            step = (T_max - T_min) / (n_replicas - 1)
            for i in range(n_replicas):
                temperatures.append(T_min + i * step)
        else:  # exponential
            # Exponential spacing
            import math
            temperatures = []
            ratio = (T_max / T_min) ** (1.0 / (n_replicas - 1))
            for i in range(n_replicas):
                temperatures.append(T_min * (ratio ** i))
        
        return temperatures
    
    def get_scaling_factors(self) -> list:
        """
        Calculate REST2 scaling factors (λ values)
        
        Returns:
            List of scaling factors for each replica
        """
        temperatures = self.get_temperature_ladder()
        T_ref = temperatures[0]  # Reference temperature (lowest)
        
        scaling_factors = []
        for T in temperatures:
            lambda_val = T_ref / T
            scaling_factors.append(lambda_val)
        
        return scaling_factors
    
    def print_summary(self) -> None:
        """Print configuration summary"""
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
        
        # Set some parameters
        config_manager.set_parameter('T_max', 450.0)
        config_manager.set_parameter('n_replicas', 12)
        config_manager.set_parameter('md_results_dir', '/path/to/md_results')
        
        # Print summary
        config_manager.print_summary()
        
        # Create and load from file
        config_manager.create_default_config("test_config.yaml")
        print("Default configuration file created: test_config.yaml")
        
    except ConfigValidationError as e:
        print(f"Configuration error: {e}")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()