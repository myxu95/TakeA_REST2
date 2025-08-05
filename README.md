# REST2 Enhanced Sampling Automation

Automated setup and execution of REST2 (Replica Exchange with Solute Tempering) enhanced sampling simulations using GROMACS and PLUMED.

## Quick Start

1. **Prepare MD results**: Complete standard MD simulation (EM → NVT → NPT → MD)
2. **Configure**: Edit `config.yaml` with your system parameters
3. **Run**: `python main.py -c config.yaml`
4. **Execute**: Submit generated scripts to run REST2 simulation

## Requirements

- Python 3.7+
- GROMACS with MPI support
- PLUMED (compiled with GROMACS)
- MDAnalysis: `pip install MDAnalysis`

## Project Structure

```
REST2_Project/
├── main.py                    # Main execution script
├── config.yaml               # Configuration file
├── modules/                   # Core modules
│   ├── config_manager.py
│   ├── structure_analyzer.py
│   ├── solute_selector.py
│   ├── replica_generator.py
│   ├── temperature_controller.py
│   └── gromacs_runner.py
├── templates/
│   └── plumed.dat            # PLUMED template (optional)
└── example/
    └── MD_results/           # Your completed MD simulation
        ├── md.tpr
        ├── topol.top
        ├── md.gro
        └── md.xtc (optional)
```

## Configuration File

Create `config.yaml` with your system parameters:

```yaml
# Basic REST2 settings
target_type: peptide          # 'peptide' or 'small_molecule'
T_min: 300.0                 # Minimum temperature (K)
T_max: 340.0                 # Maximum temperature (K)
replex: 200                  # Exchange interval (steps)
n_replicas: 8                # Number of replicas
scaling_method: linear       # 'linear' or 'exponential'

# Target selection
target_selection: chain A    # 'chain A' for peptide or 'resname LIG' for ligand
distance_range: 6.0          # Cutoff distance (Angstrom)
use_trajectory: false        # Use trajectory for dynamic selection
occupancy_threshold: 0.5     # Occupancy threshold if using trajectory

# File paths
md_results_dir: example/MD_results  # Path to completed MD simulation
plumed_dat: templates/plumed.dat    # PLUMED template (optional)
output_dir: ./rest2_simulation      # Output directory

# GROMACS execution settings
gromacs:
  gmx_mpi_command: gmx_mpi
  n_cpus: 8                  # Number of CPU cores
  n_gpus: 2                  # Number of GPUs
  script_types:              # Scripts to generate
    - slurm
    - localrun
    - test
```

## Usage

### Basic Usage
```bash
# Complete REST2 setup
python main.py -c config.yaml

# Custom output directory
python main.py -c config.yaml -o ./my_rest2_sim
```

### Advanced Usage
```bash
# Only validate configuration
python main.py -c config.yaml --validate-only

# Only generate execution scripts
python main.py -c config.yaml --scripts-only

# Verbose output
python main.py -c config.yaml --verbose
```

## Workflow Steps

The automation performs these steps:

1. **Structure Analysis**: Identifies target atoms and nearby residues
2. **Solute Selection**: Modifies topology file to mark REST2 atoms
3. **Replica Generation**: Creates replica directories and copies files
4. **Temperature Control**: Generates temperature-specific input files with PLUMED REST2
5. **Script Generation**: Creates execution scripts for your environment

## Generated Files

After successful execution:

```
rest2_simulation/
├── replica_0/
│   ├── input/
│   │   ├── input.tpr         # GROMACS run file
│   │   ├── topol.top         # Modified topology
│   │   ├── rest2.mdp         # MD parameters
│   │   ├── plumed.dat        # PLUMED with PARTIAL_TEMPERING
│   │   └── index.ndx         # Index file
│   └── output/               # Simulation outputs (created during run)
├── replica_1/
├── ...
├── run_rest2.slurm          # SLURM job script
├── run_rest2_local.sh       # Local execution script
├── test_rest2.sh            # Quick test script
└── temperature_summary.txt  # Temperature and scaling info
```

## Running Simulations

### Test Setup
```bash
cd rest2_simulation
./test_rest2.sh              # Quick validation (1000 steps)
```

### Cluster Submission
```bash
sbatch run_rest2.slurm        # Submit to SLURM scheduler
```

### Local Execution
```bash
./run_rest2_local.sh          # Run on workstation
```

## Configuration Options

### Target Types
- **peptide**: Use `target_selection: chain A` (or other chain)
- **small_molecule**: Use `target_selection: resname LIG` (or other residue name)

### Temperature Scaling
- **linear**: Even temperature spacing
- **exponential**: Exponential temperature distribution

### Resource Configuration
```yaml
gromacs:
  n_cpus: 8                   # Total CPU cores
  n_gpus: 2                   # Total GPUs
  script_types:               # Choose scripts to generate
    - slurm                   # Cluster submission
    - localrun                # Local execution
    - test                    # Quick validation
```

### Trajectory-based Selection
```yaml
use_trajectory: true          # Use MD trajectory for selection
occupancy_threshold: 0.5      # Minimum contact occupancy (0-1)
```

## Input Requirements

Your `MD_results/` directory must contain:
- `md.tpr`: Production MD run file
- `topol.top`: System topology file  
- `md.gro`: Final structure file
- `md.xtc`: Trajectory file (if `use_trajectory: true`)

## PLUMED Integration

The system automatically generates PLUMED `PARTIAL_TEMPERING` commands:

```plumed
PARTIAL_TEMPERING ...
  ATOMS=1-50,75-100           # Solute atoms (auto-generated)
  TEMP=300.0                  # Reference temperature
  LAMBDA=0.882353             # Scaling factor for this replica
  LABEL=rest2_scaling
... PARTIAL_TEMPERING
```

## Troubleshooting

### Common Issues

**Configuration validation failed**
- Check file paths in `md_results_dir`
- Verify temperature range: `T_min < T_max`
- Ensure `target_selection` matches your system

**Structure analysis failed**
- Verify structure and topology files exist
- Check `target_selection` syntax (GROMACS selection language)
- Ensure trajectory file exists if `use_trajectory: true`

**No atoms found for selection**
- Check `target_selection` string
- Verify chain/residue names in your system
- Try simpler selection like `protein` or `not water`

**Script execution fails**
- Verify GROMACS and PLUMED installation
- Check MPI configuration
- Ensure sufficient GPU/CPU resources

### File Requirements

Ensure your MD results contain standard GROMACS files:
```bash
ls example/MD_results/
# Should show: md.tpr topol.top md.gro [md.xtc]
```

## Performance Tips

- Use 1 GPU per replica for optimal performance
- Set `n_cpus` equal to `n_replicas` for standard setups
- For large systems, consider fewer replicas with smaller temperature gaps
- Test with short runs before full production simulations

## Support

For issues or questions:
1. Check configuration file syntax
2. Verify input file requirements
3. Test with `--validate-only` flag
4. Use `--verbose` for detailed error messages