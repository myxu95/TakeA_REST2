#!/usr/bin/env python3
"""
REST2 Enhanced Sampling Project - Main Script
Automated REST2 simulation setup and execution
"""

import sys
import argparse
from pathlib import Path

# Import all project modules
from modules.config_manager import ConfigManager, ConfigValidationError
from modules.structure_analyzer import StructureAnalyzer, StructureAnalysisError
from modules.solute_selector import SoluteSelector, SoluteSelectorError
from modules.replica_generator import ReplicaGenerator, ReplicaGeneratorError
from modules.temperature_controller import TemperatureController, TemperatureControllerError
from modules.gromacs_runner import GromacsRunner


def setup_argument_parser():
    """Setup command line argument parser"""
    parser = argparse.ArgumentParser(
        description="REST2 Enhanced Sampling Automation Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Complete setup with default config
  python main.py -c config.yaml

  # Setup with custom output directory
  python main.py -c config.yaml -o ./my_rest2_simulation

  # Only generate scripts (skip analysis and setup)
  python main.py -c config.yaml --scripts-only

  # Validate configuration without running
  python main.py -c config.yaml --validate-only
        """
    )
    
    parser.add_argument(
        '-c', '--config', 
        required=True,
        help='Configuration file path'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        help='Output directory (overrides config)'
    )
    
    parser.add_argument(
        '--scripts-only',
        action='store_true',
        help='Only generate execution scripts (requires existing setup)'
    )
    
    parser.add_argument(
        '--validate-only',
        action='store_true',
        help='Only validate configuration without execution'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    return parser


def validate_configuration(config_manager):
    """Validate configuration and print summary"""
    print("="*60)
    print("Configuration Validation")
    print("="*60)
    
    try:
        config_manager.validate_config()
        config_manager.print_summary()
        print("✓ Configuration validation passed")
        return True
    except ConfigValidationError as e:
        print(f"✗ Configuration validation failed:")
        print(f"  {e}")
        return False


def run_structure_analysis(config_manager, verbose=False):
    """Run structure analysis step"""
    print("\n" + "="*60)
    print("Step 1: Structure Analysis")
    print("="*60)
    
    try:
        # Get file paths from config
        md_results_dir = config_manager.get_parameter('md_results_dir')
        structure_file = str(Path(md_results_dir) / 'md.gro')
        topology_file = config_manager.get_parameter('input_tpr')
        trajectory_file = str(Path(md_results_dir) / 'md.xtc') if config_manager.get_parameter('use_trajectory') else None
        
        # Initialize analyzer
        analyzer = StructureAnalyzer(structure_file, topology_file, trajectory_file)
        
        # Run analysis
        results = analyzer.analyze_target_and_environment(
            target_selection=config_manager.get_parameter('target_selection'),
            cutoff_distance=config_manager.get_parameter('distance_range'),
            use_trajectory=config_manager.get_parameter('use_trajectory'),
            occupancy_threshold=config_manager.get_parameter('occupancy_threshold')
        )
        
        # Print results
        analyzer.print_selected_residues(results)
        
        # Get data for next step
        solute_data = analyzer.get_solute_selection_data(results)
        
        print("✓ Structure analysis completed")
        return solute_data
        
    except (StructureAnalysisError, FileNotFoundError) as e:
        print(f"✗ Structure analysis failed: {e}")
        return None


def run_solute_selection(config_manager, solute_data, verbose=False):
    """Run solute selection step"""
    print("\n" + "="*60)
    print("Step 2: Solute Selection and Topology Modification")
    print("="*60)
    
    try:
        # Initialize solute selector
        selector = SoluteSelector(solute_data)
        
        # Print modification summary
        selector.print_modification_summary()
        
        # Modify topology file
        input_topology = config_manager.get_parameter('topology')
        output_topology = config_manager.get_parameter('output_tpr', 'rest2_topol.top')
        
        selector.modify_topology_file(input_topology, output_topology)
        
        # Validate modification
        if selector.validate_topology_modification(output_topology):
            print("✓ Solute selection and topology modification completed")
            return output_topology
        else:
            print("✗ Topology modification validation failed")
            return None
            
    except SoluteSelectorError as e:
        print(f"✗ Solute selection failed: {e}")
        return None


def run_replica_generation(config_manager, modified_topology, verbose=False):
    """Run replica generation step"""
    print("\n" + "="*60)
    print("Step 3: Replica Generation")
    print("="*60)
    
    try:
        # Initialize replica generator
        generator = ReplicaGenerator(config_manager)
        
        # Print configuration summary
        generator.print_replica_summary()
        
        # Set up replica directories
        generator.setup_replica_directories()
        
        # Copy base files
        generator.copy_base_files_to_replicas(modified_topology)
        
        # Create replica info files
        generator.create_replica_info_files()
        
        # Validate setup
        if generator.validate_replica_setup():
            replica_data = generator.get_replica_data()
            print("✓ Replica generation completed")
            return replica_data
        else:
            print("✗ Replica setup validation failed")
            return None
            
    except ReplicaGeneratorError as e:
        print(f"✗ Replica generation failed: {e}")
        return None


def run_temperature_control(config_manager, replica_data, solute_data, verbose=False):
    """Run temperature control step"""
    print("\n" + "="*60)
    print("Step 4: Temperature Control and File Preparation")
    print("="*60)
    
    try:
        # Initialize temperature controller
        controller = TemperatureController(config_manager, replica_data, solute_data)
        
        # Print configuration summary
        controller.print_temperature_summary()
        
        # Generate scaled topology files
        base_topology = config_manager.get_parameter('output_tpr', 'rest2_topol.top')
        controller.generate_scaled_topology_files(base_topology)
        
        # Generate MDP files
        controller.generate_mdp_files()
        
        # Prepare additional input files
        controller.prepare_additional_input_files()
        
        # Create temperature summary
        controller.create_temperature_summary()
        
        # Validate setup
        if controller.validate_temperature_setup():
            print("✓ Temperature control and file preparation completed")
            return True
        else:
            print("✗ Temperature setup validation failed")
            return False
            
    except TemperatureControllerError as e:
        print(f"✗ Temperature control failed: {e}")
        return False


def run_script_generation(config_manager, replica_data, verbose=False):
    """Run script generation step"""
    print("\n" + "="*60)
    print("Step 5: Execution Script Generation")
    print("="*60)
    
    try:
        # Initialize GROMACS runner
        runner = GromacsRunner(config_manager, replica_data)
        
        # Create scripts based on configuration
        created_scripts = runner.create_configured_scripts()
        
        # Print summary
        runner.print_summary(created_scripts)
        
        print("✓ Execution script generation completed")
        return created_scripts
        
    except Exception as e:
        print(f"✗ Script generation failed: {e}")
        return None


def main():
    """Main execution function"""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    print("="*60)
    print("REST2 Enhanced Sampling Automation")
    print("="*60)
    print(f"Configuration file: {args.config}")
    if args.output_dir:
        print(f"Output directory: {args.output_dir}")
    print()
    
    try:
        # Load configuration
        config_manager = ConfigManager(args.config)
        
        # Override output directory if specified
        if args.output_dir:
            config_manager.set_parameter('output_dir', args.output_dir)
        
        # Validate configuration
        if not validate_configuration(config_manager):
            sys.exit(1)
        
        # If only validation requested, exit here
        if args.validate_only:
            print("\n✓ Configuration validation completed successfully")
            sys.exit(0)
        
        # If only script generation requested
        if args.scripts_only:
            print("\nGenerating execution scripts only...")
            
            # Create minimal replica data for script generation
            replica_data = {
                'n_replicas': config_manager.get_parameter('n_replicas'),
                'base_output_dir': config_manager.get_parameter('output_dir')
            }
            
            created_scripts = run_script_generation(config_manager, replica_data, args.verbose)
            
            if created_scripts:
                print("\n✓ Script generation completed successfully")
                sys.exit(0)
            else:
                print("\n✗ Script generation failed")
                sys.exit(1)
        
        # Full workflow execution
        print("Starting full REST2 setup workflow...")
        
        # Step 1: Structure Analysis
        solute_data = run_structure_analysis(config_manager, args.verbose)
        if not solute_data:
            sys.exit(1)
        
        # Step 2: Solute Selection
        modified_topology = run_solute_selection(config_manager, solute_data, args.verbose)
        if not modified_topology:
            sys.exit(1)
        
        # Step 3: Replica Generation
        replica_data = run_replica_generation(config_manager, modified_topology, args.verbose)
        if not replica_data:
            sys.exit(1)
        
        # Step 4: Temperature Control
        if not run_temperature_control(config_manager, replica_data, solute_data, args.verbose):
            sys.exit(1)
        
        # Step 5: Script Generation
        created_scripts = run_script_generation(config_manager, replica_data, args.verbose)
        if not created_scripts:
            sys.exit(1)
        
        # Final summary
        print("\n" + "="*60)
        print("REST2 Setup Completed Successfully!")
        print("="*60)
        print("All components have been generated and validated.")
        print(f"Output directory: {config_manager.get_parameter('output_dir')}")
        
        if created_scripts:
            print("\nNext steps:")
            for script_type, script_path in created_scripts.items():
                script_name = Path(script_path).name
                if script_type == 'test':
                    print(f"  1. Test setup: ./{script_name}")
                elif script_type == 'slurm':
                    print(f"  2. Submit job: sbatch {script_name}")
                elif script_type == 'localrun':
                    print(f"  2. Run locally: ./{script_name}")
        
        print("\n✓ REST2 Enhanced Sampling setup completed successfully!")
        
    except KeyboardInterrupt:
        print("\n\nSetup interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()