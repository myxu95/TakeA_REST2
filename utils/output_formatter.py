#!/usr/bin/env python3
"""
REST2 Output Formatter Module
Unified output formatting for REST2 enhanced sampling simulations
"""

import sys
from typing import Dict, List, Any, Optional, Union
from datetime import datetime


class OutputFormatter:
    """
    Unified output formatter for REST2 simulations
    Handles all output formatting across different modules
    """
    
    # Default formatting styles
    STYLES = {
        'header': '=' * 60,
        'subheader': '-' * 50,
        'separator': '-' * 40,
        'success': '✓',
        'error': '✗',
        'warning': '⚠',
        'info': 'ℹ'
    }
    
    @staticmethod
    def print_header(title: str, style: str = 'header') -> None:
        """
        Print formatted header
        
        Args:
            title: Header title
            style: Header style
        """
        print(f"\n{OutputFormatter.STYLES[style]}")
        print(f"{title}")
        print(f"{OutputFormatter.STYLES[style]}")
    
    @staticmethod
    def print_subheader(title: str, style: str = 'subheader') -> None:
        """
        Print formatted subheader
        
        Args:
            title: Subheader title
            style: Subheader style
        """
        print(f"\n{title}")
        print(f"{OutputFormatter.STYLES[style]}")
    
    @staticmethod
    def print_section(title: str, content: Union[str, List[str]], 
                     style: str = 'separator') -> None:
        """
        Print formatted section
        
        Args:
            title: Section title
            content: Section content
            style: Section style
        """
        print(f"\n{title}")
        print(f"{OutputFormatter.STYLES[style]}")
        
        if isinstance(content, str):
            print(content)
        elif isinstance(content, list):
            for item in content:
                print(f"  {item}")
    
    @staticmethod
    def print_summary(data: Dict[str, Any], title: str = "Summary") -> None:
        """
        Print formatted summary
        
        Args:
            data: Summary data dictionary
            title: Summary title
        """
        OutputFormatter.print_subheader(title)
        
        for key, value in data.items():
            if isinstance(value, (int, float)):
                print(f"{key:<25}: {value}")
            elif isinstance(value, list):
                print(f"{key:<25}: {len(value)} items")
                for i, item in enumerate(value[:5]):  # Show first 5 items
                    print(f"  {'':<25}  {i}: {item}")
                if len(value) > 5:
                    print(f"  {'':<25}  ... and {len(value) - 5} more")
            elif isinstance(value, dict):
                print(f"{key:<25}: {len(value)} keys")
                for sub_key, sub_value in list(value.items())[:3]:  # Show first 3
                    print(f"  {'':<25}  {sub_key}: {sub_value}")
                if len(value) > 3:
                    print(f"  {'':<25}  ... and {len(value) - 3} more")
            else:
                print(f"{key:<25}: {value}")
    
    @staticmethod
    def print_table(headers: List[str], rows: List[List[Any]], 
                   title: str = "Table") -> None:
        """
        Print formatted table
        
        Args:
            headers: Table headers
            rows: Table rows
            title: Table title
        """
        if not rows:
            return
        
        OutputFormatter.print_subheader(title)
        
        # Calculate column widths
        col_widths = []
        for i, header in enumerate(headers):
            max_width = len(str(header))
            for row in rows:
                if i < len(row):
                    max_width = max(max_width, len(str(row[i])))
            col_widths.append(max_width)
        
        # Print header
        header_str = "  ".join(f"{header:<{width}}" for header, width in zip(headers, col_widths))
        print(header_str)
        print("-" * len(header_str))
        
        # Print rows
        for row in rows:
            row_str = "  ".join(f"{str(cell):<{width}}" for cell, width in zip(row, col_widths))
            print(row_str)
    
    @staticmethod
    def print_status(message: str, status: str = 'info', indent: int = 0) -> None:
        """
        Print status message
        
        Args:
            message: Status message
            status: Status type ('success', 'error', 'warning', 'info')
            indent: Indentation level
        """
        icon = OutputFormatter.STYLES.get(status, '')
        indent_str = "  " * indent
        print(f"{indent_str}{icon} {message}")
    
    @staticmethod
    def print_progress(current: int, total: int, description: str = "Progress") -> None:
        """
        Print progress bar
        
        Args:
            current: Current progress
            total: Total steps
            description: Progress description
        """
        if total == 0:
            return
        
        percentage = (current / total) * 100
        bar_length = 30
        filled_length = int(bar_length * current // total)
        bar = '█' * filled_length + '░' * (bar_length - filled_length)
        
        print(f"\r{description}: |{bar}| {percentage:.1f}% ({current}/{total})", end='')
        sys.stdout.flush()
        
        if current == total:
            print()  # New line when complete
    
    @staticmethod
    def print_validation_summary(errors: List[str], context: str = "Validation") -> None:
        """
        Print validation summary
        
        Args:
            errors: List of validation errors
            context: Validation context
        """
        if not errors:
            OutputFormatter.print_status(f"{context} passed", 'success')
        else:
            OutputFormatter.print_status(f"{context} failed", 'error')
            for error in errors:
                OutputFormatter.print_status(f"  - {error}", 'error', indent=1)
    
    @staticmethod
    def print_configuration_summary(config: Dict[str, Any]) -> None:
        """
        Print configuration summary
        
        Args:
            config: Configuration dictionary
        """
        OutputFormatter.print_header("REST2 Configuration Summary")
        
        # Basic settings
        basic_settings = {
            'Target type': config.get('target_type', 'Not specified'),
            'Target selection': config.get('target_selection', 'Not specified'),
            'Temperature range': f"{config.get('T_min', 0):.1f} - {config.get('T_max', 0):.1f} K",
            'Number of replicas': config.get('n_replicas', 0),
            'Exchange interval': f"{config.get('replex', 0)} steps",
            'Scaling method': config.get('scaling_method', 'Not specified'),
            'Distance cutoff': f"{config.get('distance_range', 0)} Angstrom"
        }
        
        OutputFormatter.print_summary(basic_settings, "Basic Settings")
        
        # Selection settings
        selection_settings = {}
        if config.get('use_trajectory', False):
            selection_settings['Selection method'] = 'Trajectory-based'
            selection_settings['Occupancy threshold'] = config.get('occupancy_threshold', 0.5)
        else:
            selection_settings['Selection method'] = 'Static structure'
        
        if selection_settings:
            OutputFormatter.print_summary(selection_settings, "Selection Settings")
        
        # File settings
        file_settings = {
            'TPR file': config.get('input_tpr', 'Not specified'),
            'Topology': config.get('topology', 'Not specified'),
            'PLUMED': config.get('plumed_dat', 'Not specified'),
            'Output directory': config.get('output_dir', 'Not specified')
        }
        
        OutputFormatter.print_summary(file_settings, "File Settings")
    
    @staticmethod
    def print_temperature_summary(temperatures: List[float], scaling_factors: List[float],
                                method: str = 'linear') -> None:
        """
        Print temperature ladder summary
        
        Args:
            temperatures: List of temperatures
            scaling_factors: List of scaling factors
            method: Scaling method
        """
        OutputFormatter.print_subheader(f"Temperature Ladder ({method} scaling)")
        
        # Create table data
        headers = ['Replica', 'Temperature (K)', 'Scaling Factor (λ)']
        rows = []
        
        for i, (T, lambda_val) in enumerate(zip(temperatures, scaling_factors)):
            rows.append([i, f"{T:.1f}", f"{lambda_val:.6f}"])
        
        OutputFormatter.print_table(headers, rows)
        
        # Print summary info
        summary_data = {
            'Temperature range': f"{min(temperatures):.1f} - {max(temperatures):.1f} K",
            'Reference temperature': f"{min(temperatures):.1f} K",
            'Number of replicas': len(temperatures)
        }
        
        OutputFormatter.print_summary(summary_data, "Temperature Summary")
    
    @staticmethod
    def print_replica_summary(replica_data: Dict[str, Any]) -> None:
        """
        Print replica setup summary
        
        Args:
            replica_data: Replica data dictionary
        """
        OutputFormatter.print_header("Replica Generation Summary")
        
        # Basic replica info
        basic_info = {
            'Number of replicas': replica_data.get('n_replicas', 0),
            'Output directory': replica_data.get('base_output_dir', 'Not specified'),
            'Scaling method': replica_data.get('scaling_method', 'Not specified')
        }
        
        OutputFormatter.print_summary(basic_info, "Replica Configuration")
        
        # Directory structure
        n_replicas = replica_data.get('n_replicas', 0)
        structure_info = []
        for i in range(n_replicas):
            structure_info.append(f"replica_{i}/")
            structure_info.append(f"  ├── input/     (TPR, topology, PLUMED)")
            structure_info.append(f"  ├── output/    (simulation results)")
            structure_info.append(f"  └── replica_info.txt")
        
        OutputFormatter.print_section("Directory Structure", structure_info)
    
    @staticmethod
    def print_execution_summary(scripts: Dict[str, str], config: Dict[str, Any]) -> None:
        """
        Print execution script summary
        
        Args:
            scripts: Dictionary of created scripts
            config: Configuration dictionary
        """
        OutputFormatter.print_header("Execution Script Generation Summary")
        
        # Script summary
        script_info = {
            'Total scripts created': len(scripts),
            'Script types': list(scripts.keys())
        }
        
        OutputFormatter.print_summary(script_info, "Script Generation")
        
        # Next steps
        next_steps = []
        for script_type, script_path in scripts.items():
            script_name = script_path.split('/')[-1]
            if script_type == 'test':
                next_steps.append(f"1. Test setup: ./{script_name}")
            elif script_type == 'slurm':
                next_steps.append(f"2. Submit job: sbatch {script_name}")
            elif script_type == 'localrun':
                next_steps.append(f"2. Run locally: ./{script_name}")
        
        if next_steps:
            OutputFormatter.print_section("Next Steps", next_steps)
    
    @staticmethod
    def print_complete_summary(config: Dict[str, Any], replica_data: Optional[Dict[str, Any]] = None,
                             solute_data: Optional[Dict[str, Any]] = None,
                             scripts: Optional[Dict[str, str]] = None) -> None:
        """
        Print complete REST2 setup summary
        
        Args:
            config: Configuration dictionary
            replica_data: Replica data dictionary (optional)
            solute_data: Solute data dictionary (optional)
            scripts: Script data dictionary (optional)
        """
        OutputFormatter.print_header("REST2 Setup Completed Successfully!")
        
        # Configuration summary
        OutputFormatter.print_configuration_summary(config)
        
        # Replica summary if available
        if replica_data:
            OutputFormatter.print_replica_summary(replica_data)
        
        # Solute summary if available
        if solute_data:
            solute_summary = {
                'Target residues': len(solute_data.get('target_residues', [])),
                'Nearby residues': len(solute_data.get('nearby_residues', [])),
                'Total solute residues': solute_data.get('total_residues', 0),
                'Target atoms': solute_data.get('target_atoms', 0),
                'Solute atoms': solute_data.get('solute_atoms', 0)
            }
            OutputFormatter.print_summary(solute_summary, "Solute Selection")
        
        # Script summary if available
        if scripts:
            OutputFormatter.print_execution_summary(scripts, config)
        
        # Final message
        OutputFormatter.print_status("All components have been generated and validated.", 'success')
        OutputFormatter.print_status(f"Output directory: {config.get('output_dir', 'Not specified')}", 'info')


def main():
    """Test output formatter"""
    try:
        print("Testing Output Formatter")
        print("=" * 40)
        
        # Test basic formatting
        OutputFormatter.print_header("Test Header")
        OutputFormatter.print_subheader("Test Subheader")
        OutputFormatter.print_section("Test Section", ["Item 1", "Item 2", "Item 3"])
        
        # Test status messages
        OutputFormatter.print_status("Success message", 'success')
        OutputFormatter.print_status("Error message", 'error')
        OutputFormatter.print_status("Warning message", 'warning')
        OutputFormatter.print_status("Info message", 'info')
        
        # Test table
        headers = ['Name', 'Age', 'City']
        rows = [
            ['Alice', 25, 'New York'],
            ['Bob', 30, 'Los Angeles'],
            ['Charlie', 35, 'Chicago']
        ]
        OutputFormatter.print_table(headers, rows, "Test Table")
        
        # Test summary
        test_data = {
            'String': 'Test value',
            'Number': 42,
            'List': ['item1', 'item2', 'item3', 'item4', 'item5', 'item6'],
            'Dict': {'key1': 'value1', 'key2': 'value2', 'key3': 'value3'}
        }
        OutputFormatter.print_summary(test_data, "Test Summary")
        
        # Test configuration summary
        test_config = {
            'target_type': 'peptide',
            'target_selection': 'chain A',
            'T_min': 300.0,
            'T_max': 340.0,
            'n_replicas': 8,
            'replex': 200,
            'scaling_method': 'linear',
            'distance_range': 6.0,
            'use_trajectory': False,
            'input_tpr': 'example/md.tpr',
            'topology': 'example/topol.top',
            'plumed_dat': 'templates/plumed.dat',
            'output_dir': './rest2_simulation'
        }
        OutputFormatter.print_configuration_summary(test_config)
        
        # Test temperature summary
        temperatures = [300.0, 305.7, 311.4, 317.1, 322.9, 328.6, 334.3, 340.0]
        scaling_factors = [1.0, 0.981, 0.963, 0.946, 0.929, 0.913, 0.897, 0.882]
        OutputFormatter.print_temperature_summary(temperatures, scaling_factors, 'linear')
        
        print("\n✓ Output formatter tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 