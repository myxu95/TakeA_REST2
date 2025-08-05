#!/usr/bin/env python3
"""
REST2 GROMACS Runner Module
Generates execution scripts for REST2 enhanced sampling simulations
"""

import os
from pathlib import Path
from typing import Dict, Any


class GromacsRunner:
    """
    GROMACS runner for REST2 simulations
    Generates execution scripts for different environments
    """
    
    def __init__(self, config_manager, replica_data: Dict[str, Any]):
        """
        Initialize GROMACS runner
        
        Args:
            config_manager: ConfigManager instance
            replica_data: Data from ReplicaGenerator
        """
        self.config = config_manager
        self.n_replicas = replica_data['n_replicas']
        self.base_output_dir = Path(replica_data['base_output_dir'])
        
        # GROMACS settings
        self.gmx_mpi_command = config_manager.get_parameter('gromacs.gmx_mpi_command', 'gmx_mpi')
        self.replex = config_manager.get_parameter('replex')
        self.plumed_dat = config_manager.get_parameter('plumed_dat', '')
        
        # Resource settings
        self.n_cpus = config_manager.get_parameter('gromacs.n_cpus', self.n_replicas)
        self.n_gpus = config_manager.get_parameter('gromacs.n_gpus', self.n_replicas)
        
        # Script generation settings
        self.script_types = config_manager.get_parameter('gromacs.script_types', ['slurm', 'localrun', 'test'])
        
        # Simulation settings
        production_time_ns = config_manager.get_parameter('simulation.production_time', 100.0)
        self.nsteps = int(production_time_ns * 1000 / 0.002)  # Convert ns to steps
        
        # Multidir string
        self.multidir_string = " ".join([f"replica_{i}" for i in range(self.n_replicas)])
    
    def create_slurm_script(self, script_name: str = "run_rest2.slurm") -> str:
        """Create SLURM job submission script"""
        script_path = self.base_output_dir / script_name
        
        # PLUMED option
        plumed_option = f" -plumed input/{Path(self.plumed_dat).name}" if self.plumed_dat else ""
        
        script_content = f"""#!/bin/bash
#SBATCH --job-name=REST2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={self.n_cpus}
#SBATCH --gres=gpu:{self.n_gpus}
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G

# Load modules (adjust for your cluster)
# module load gromacs/2023
# module load plumed/2.8

export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR

# Set GPU IDs (use specified number of GPUs)
if [ {self.n_gpus} -gt 1 ]; then
    gpu_ids=$(seq -s "," 0 $(({self.n_gpus}-1)))
else
    gpu_ids="0"
fi

echo "Starting REST2 simulation:"
echo "  Replicas: {self.n_replicas}"
echo "  CPUs: {self.n_cpus}"
echo "  GPUs: {self.n_gpus} (IDs: $gpu_ids)"
echo "  Exchange interval: {self.replex} steps"

mpirun -np {self.n_replicas} --oversubscribe {self.gmx_mpi_command} mdrun \\
    -v -deffnm output/rest2 \\
    -s input/input.tpr \\
    -multidir {self.multidir_string} \\
    -replex {self.replex} \\
    -nb gpu -bonded gpu -pme gpu -update gpu \\
    -npme -1 -hrex{plumed_option} \\
    -nsteps {self.nsteps} \\
    -gpu_id "$gpu_ids"

echo "REST2 simulation completed at $(date)"
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def create_localrun_script(self, script_name: str = "run_rest2_local.sh") -> str:
        """Create local run script"""
        script_path = self.base_output_dir / script_name
        
        # PLUMED option
        plumed_option = f" -plumed input/{Path(self.plumed_dat).name}" if self.plumed_dat else ""
        
        script_content = f"""#!/bin/bash

echo "REST2 Local Run - {self.n_replicas} replicas"
echo "CPUs: {self.n_cpus}, GPUs: {self.n_gpus}"
echo "Start time: $(date)"

# Check GROMACS
if ! command -v {self.gmx_mpi_command} &> /dev/null; then
    echo "Error: {self.gmx_mpi_command} not found"
    exit 1
fi

# Set GPU configuration
if [ {self.n_gpus} -gt 0 ] && command -v nvidia-smi &> /dev/null; then
    available_gpus=$(nvidia-smi --list-gpus | wc -l)
    if [ $available_gpus -ge {self.n_gpus} ]; then
        if [ {self.n_gpus} -gt 1 ]; then
            gpu_ids=$(seq -s "," 0 $(({self.n_gpus}-1)))
        else
            gpu_ids="0"
        fi
        echo "Using GPUs: $gpu_ids"
        USE_GPU=true
    else
        echo "Warning: Only $available_gpus GPUs available, requested {self.n_gpus}"
        gpu_ids="0"
        USE_GPU=true
    fi
else
    echo "No GPU usage configured or no GPUs detected"
    USE_GPU=false
fi

# Validate replica directories
for i in $(seq 0 {self.n_replicas-1}); do
    if [ ! -f "replica_$i/input/input.tpr" ]; then
        echo "Error: replica_$i/input/input.tpr not found"
        exit 1
    fi
done

echo "Starting simulation..."

if [ "$USE_GPU" = true ]; then
    mpirun -np {self.n_replicas} --oversubscribe {self.gmx_mpi_command} mdrun \\
        -v -deffnm output/rest2 \\
        -s input/input.tpr \\
        -multidir {self.multidir_string} \\
        -replex {self.replex} \\
        -nb gpu -bonded gpu -pme gpu -update gpu \\
        -npme -1 -hrex{plumed_option} \\
        -nsteps {self.nsteps} \\
        -gpu_id "$gpu_ids"
else
    mpirun -np {self.n_replicas} --oversubscribe {self.gmx_mpi_command} mdrun \\
        -v -deffnm output/rest2 \\
        -s input/input.tpr \\
        -multidir {self.multidir_string} \\
        -replex {self.replex} \\
        -hrex{plumed_option} \\
        -nsteps {self.nsteps}
fi

if [ $? -eq 0 ]; then
    echo "Simulation completed successfully at $(date)"
else
    echo "Simulation failed"
    exit 1
fi
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def create_test_script(self, script_name: str = "test_rest2.sh") -> str:
        """Create quick test script"""
        script_path = self.base_output_dir / script_name
        
        plumed_option = f" -plumed input/{Path(self.plumed_dat).name}" if self.plumed_dat else ""
        
        script_content = f"""#!/bin/bash

echo "REST2 Test Run (1000 steps)"

# Quick validation
for i in $(seq 0 {self.n_replicas-1}); do
    if [ ! -f "replica_$i/input/input.tpr" ]; then
        echo "Error: replica_$i/input/input.tpr not found"
        exit 1
    fi
done

mpirun -np {self.n_replicas} --oversubscribe {self.gmx_mpi_command} mdrun \\
    -v -deffnm output/test \\
    -s input/input.tpr \\
    -multidir {self.multidir_string} \\
    -replex {self.replex} \\
    -hrex{plumed_option} \\
    -nsteps 1000

if [ $? -eq 0 ]; then
    echo "Test completed successfully!"
else
    echo "Test failed"
    exit 1
fi
"""
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        
        return str(script_path)
    
    def create_scripts(self, script_types: list = None) -> Dict[str, str]:
        """
        Create selected execution scripts
        
        Args:
            script_types: List of script types to create ('slurm', 'localrun', 'test')
                         If None, creates all scripts
        
        Returns:
            Dictionary with script types and their file paths
        """
        if script_types is None:
            script_types = ['slurm', 'localrun', 'test']
        
        created_scripts = {}
        
        for script_type in script_types:
            if script_type == 'slurm':
                created_scripts['slurm'] = self.create_slurm_script()
            elif script_type == 'localrun':
                created_scripts['localrun'] = self.create_localrun_script()
            elif script_type == 'test':
                created_scripts['test'] = self.create_test_script()
            else:
                print(f"Warning: Unknown script type '{script_type}' ignored")
        
        return created_scripts
    
    def create_configured_scripts(self) -> Dict[str, str]:
        """
        Create scripts based on configuration settings
        
        Returns:
            Dictionary with script types and their file paths
        """
        return self.create_scripts(self.script_types)
    
    def create_all_scripts(self) -> Dict[str, str]:
        """Create all execution scripts (backward compatibility)"""
        return self.create_scripts(['slurm', 'localrun', 'test'])
    
    def print_summary(self, created_scripts: Dict[str, str] = None) -> None:
        """Print script generation summary"""
        print("\n" + "="*50)
        print("REST2 Execution Scripts Generated")
        print("="*50)
        print(f"Replicas: {self.n_replicas}")
        print(f"CPUs: {self.n_cpus}")
        print(f"GPUs: {self.n_gpus}")
        print(f"Exchange interval: {self.replex} steps")
        print(f"Steps: {self.nsteps}")
        print(f"GROMACS: {self.gmx_mpi_command}")
        if self.plumed_dat:
            print(f"PLUMED: {Path(self.plumed_dat).name}")
        
        if created_scripts:
            print("\nFiles created:")
            for script_type, script_path in created_scripts.items():
                script_name = Path(script_path).name
                if script_type == 'slurm':
                    print(f"  - {script_name} (cluster submission)")
                elif script_type == 'localrun':
                    print(f"  - {script_name} (workstation)")
                elif script_type == 'test':
                    print(f"  - {script_name} (validation)")
        
        print("="*50 + "\n")


def main():
    """Test GROMACS runner"""
    try:
        class MockConfig:
            def __init__(self):
                self.params = {
                    'gromacs.gmx_mpi_command': 'gmx_mpi',
                    'gromacs.n_cpus': 8,
                    'gromacs.n_gpus': 2,
                    'replex': 200,
                    'plumed_dat': 'templates/plumed.dat',
                    'simulation.production_time': 100.0
                }
            
            def get_parameter(self, key, default=None):
                return self.params.get(key, default)
        
        mock_replica_data = {
            'n_replicas': 4,
            'base_output_dir': './test_rest2'
        }
        
        runner = GromacsRunner(MockConfig(), mock_replica_data)
        
        # Example 1: Create all scripts
        all_scripts = runner.create_all_scripts()
        runner.print_summary(all_scripts)
        
        # Example 2: Create only SLURM script
        print("Creating only SLURM script:")
        slurm_only = runner.create_scripts(['slurm'])
        runner.print_summary(slurm_only)
        
        # Example 3: Create SLURM and test scripts
        print("Creating SLURM and test scripts:")
        selected_scripts = runner.create_scripts(['slurm', 'test'])
        runner.print_summary(selected_scripts)
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()