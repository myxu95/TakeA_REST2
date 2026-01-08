#!/usr/bin/env python3
"""
REST2 TPR Generator Module
Generates TPR files for each replica using gmx_mpi grompp to ensure version compatibility
"""

import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, Optional


class TPRGeneratorError(Exception):
    """TPR generator specific errors"""
    pass


class TPRGenerator:
    """
    TPR file generator for REST2 replicas
    Uses gmx_mpi grompp to regenerate TPR files with correct GROMACS version
    """

    def __init__(self, config_manager, replica_data: Dict[str, Any], structure_file: str):
        """
        Initialize TPR generator

        Args:
            config_manager: Configuration manager instance
            replica_data: Replica data dictionary containing:
                - n_replicas: Number of replicas
                - replicas: List of replica dictionaries with 'index' and 'input_dir'
            structure_file: Path to equilibrated structure file (npt.gro)

        Raises:
            FileNotFoundError: If structure file doesn't exist
        """
        self.config = config_manager
        self.replica_data = replica_data
        self.structure_file = Path(structure_file)

        # Get GROMACS command from config or use default
        gromacs_config = config_manager.get_parameter('gromacs', {})
        self.gmx_command = gromacs_config.get('gmx_mpi_command', 'gmx_mpi')

        # Validate structure file
        if not self.structure_file.exists():
            raise FileNotFoundError(f"Structure file not found: {self.structure_file}")

    def generate_tpr_files(self) -> bool:
        """
        Generate TPR files for all replicas

        Returns:
            bool: True if all TPR files generated successfully

        Raises:
            TPRGeneratorError: If any TPR generation fails
        """
        n_replicas = self.replica_data['n_replicas']
        replicas = self.replica_data['replicas']

        print(f"\nGenerating TPR files for {n_replicas} replicas...")

        for replica in replicas:
            replica_index = replica['index']
            print(f"  Replica {replica_index}: Generating TPR file...")

            try:
                self._generate_single_tpr(replica)
                print(f"  Replica {replica_index}: ✓ TPR file created")
            except Exception as e:
                raise TPRGeneratorError(
                    f"Failed to generate TPR for replica {replica_index}: {e}"
                )

        return True

    def _generate_single_tpr(self, replica: Dict[str, Any]) -> None:
        """
        Generate TPR file for a single replica

        Args:
            replica: Replica dictionary with 'index' and 'input_dir'

        Raises:
            FileNotFoundError: If required input files are missing
            RuntimeError: If grompp command fails
        """
        input_dir = Path(replica['input_dir'])

        # Check required files
        mdp_file = input_dir / "rest2.mdp"
        top_file = input_dir / "topol.top"
        output_tpr = input_dir / "input.tpr"

        if not mdp_file.exists():
            raise FileNotFoundError(f"MDP file not found: {mdp_file}")
        if not top_file.exists():
            raise FileNotFoundError(f"Topology file not found: {top_file}")

        # Prepare grompp command
        cmd = [
            self.gmx_command, 'grompp',
            '-f', str(mdp_file),
            '-c', str(self.structure_file),
            '-p', str(top_file),
            '-o', str(output_tpr),
            '-maxwarn', '3'
        ]

        # Run grompp
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            # Check return code
            if result.returncode != 0:
                error_msg = f"grompp failed for replica {replica['index']}:\n"
                error_msg += f"Command: {' '.join(cmd)}\n"
                error_msg += f"Return code: {result.returncode}\n"
                error_msg += f"STDERR:\n{result.stderr}\n"
                error_msg += f"STDOUT:\n{result.stdout}"
                raise RuntimeError(error_msg)

            # Verify TPR file was created
            if not output_tpr.exists():
                raise FileNotFoundError(
                    f"TPR file was not created despite successful grompp: {output_tpr}"
                )

        except FileNotFoundError as e:
            # Re-raise file not found errors
            raise
        except Exception as e:
            # Catch other exceptions (e.g., command not found)
            raise RuntimeError(
                f"Error running grompp command: {e}\n"
                f"Command: {' '.join(cmd)}\n"
                f"Make sure {self.gmx_command} is installed and in PATH"
            )

    def print_generation_summary(self) -> None:
        """Print TPR generation configuration summary"""
        print("\n" + "="*60)
        print("TPR File Generation Summary")
        print("="*60)
        print(f"GROMACS command      : {self.gmx_command}")
        print(f"Structure file       : {self.structure_file}")
        print(f"Number of replicas   : {self.replica_data['n_replicas']}")
        print("="*60)

    def validate_tpr_files(self) -> bool:
        """
        Validate that all TPR files have been correctly generated

        Returns:
            bool: True if all TPR files are valid
        """
        replicas = self.replica_data['replicas']
        all_valid = True

        for replica in replicas:
            input_dir = Path(replica['input_dir'])
            tpr_file = input_dir / "input.tpr"
            replica_index = replica['index']

            # Check existence
            if not tpr_file.exists():
                print(f"✗ TPR file not found for replica {replica_index}: {tpr_file}")
                all_valid = False
                continue

            # Check file size (TPR files are typically > 1 MB)
            file_size = tpr_file.stat().st_size
            if file_size < 1000:  # Less than 1 KB is definitely wrong
                print(f"✗ TPR file too small for replica {replica_index}: {file_size} bytes")
                all_valid = False
                continue

        if all_valid:
            print("✓ All TPR files validated successfully")

        return all_valid


def main():
    """
    Test TPR generator functionality
    """
    print("TPR Generator Module Test")
    print("This module is meant to be imported, not run directly")
    print("See main.py for usage example")


if __name__ == '__main__':
    main()
