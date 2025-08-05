# REST2 Enhanced Sampling Automation

Automated setup and execution of REST2 (Replica Exchange with Solute Tempering) enhanced sampling simulations using GROMACS and PLUMED.

## Quick Start

1. **Prepare MD results**: Complete standard MD simulation (EM → NVT → NPT → MD)
2. **Configure**: Edit `configs/config_simple.yaml` with your system parameters
3. **Run**: `python main.py -c configs/config_simple.yaml`
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
├── configs/                   # Configuration files
│   ├── config_simple.yaml    # Simple configuration template
│   ├── config_template.yaml  # Full configuration template
│   └── example_config_chain_c.yaml  # Example configuration
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

## Running the Example

The project includes a complete example to demonstrate the REST2 workflow. Follow these steps to run the example:

### Step 1: Check Example Files
```bash
# Verify example files are present
ls example/
# Should show: md.gro, md.tpr, topol.top, example.xtc, charmm36-jul2021.ff/, mdp/
```

### Step 2: Install Dependencies
```bash
# Run the installation script
./install.sh

# Or manually install Python dependencies
pip install -r requirements.txt
```

### Step 3: Configure for Example
The example uses chain C selection. The configuration is already set up in `configs/example_config_chain_c.yaml`:

```yaml
# Input files for example
topology: "example/topol.top"
input_tpr: "example/md.tpr"
md_results_dir: "example"
output_dir: "rest2_simulation"

# Target selection for chain C
target_selection: "chainid C"
distance_range: 5.0

# REST2 parameters
n_replicas: 8
temperature:
  min: 300.0
  max: 400.0
temperature_method: "linear"
replex: 1000

# GROMACS settings
gromacs:
  gmx_mpi_command: "gmx_mpi"
  n_cpus: 4
  n_gpus: 1
  script_types: ["slurm", "localrun", "test"]
```

### Step 4: Run the Example
```bash
# Run complete REST2 setup
python main.py -c configs/example_config_chain_c.yaml

# Or run with verbose output
python main.py -c configs/example_config_chain_c.yaml --verbose
```

### Step 5: Verify Output
After successful execution, check the generated files:

```bash
# Check main output directory
ls rest2_simulation/

# Check replica directories
ls rest2_simulation/replica_0/input/

# Check generated scripts
ls rest2_simulation/*.sh
ls rest2_simulation/*.slurm
```

### Step 6: Test the Setup
```bash
# Navigate to output directory
cd rest2_simulation

# Run quick test (1000 steps)
./test_rest2.sh

# Check test results
ls replica_*/output/
```

### Step 7: Run Full Simulation (Optional)
```bash
# For cluster submission
sbatch run_rest2.slurm

# For local execution
./run_rest2_local.sh
```

## Example Workflow Details

The example demonstrates these workflow steps:

1. **Structure Analysis**: Analyzes `example/md.gro` and `example/md.tpr`
   - Extracts PDB with chain information using `gmx trjconv`
   - Identifies chain C atoms as target
   - Finds nearby residues within 5.0 Å

2. **Topology Merge**: Merges topology files using `gmx grompp -pp`
   - Processes `example/topol.top` with `example/md.gro`
   - Creates `processed.top` with merged topology information

3. **Solute Selection**: Modifies topology for REST2 scaling
   - Marks all solute atoms (not just alpha carbons)
   - Creates `rest2_topol.top` with REST2 modifications

4. **Replica Generation**: Creates 8 replica directories
   - `replica_0/` through `replica_7/`
   - Each with input files and output directories

5. **Temperature Control**: Generates temperature-specific files
   - Temperature range: 300K to 400K (linear spacing)
   - PLUMED `PARTIAL_TEMPERING` commands for each replica
   - MDP files with appropriate parameters

6. **Script Generation**: Creates execution scripts
   - `test_rest2.sh`: Quick validation script
   - `run_rest2_local.sh`: Local execution script
   - `run_rest2.slurm`: SLURM cluster script

## Configuration Files

### Simple Configuration (`configs/config_simple.yaml`)
Use this for basic REST2 setup:

```yaml
# Input files
topology: "example/topol.top"
input_tpr: "example/md.tpr"
md_results_dir: "example"
output_dir: "rest2_simulation"

# Target selection for REST2 scaling
target_selection: "chainid A"
distance_range: 5.0

# REST2 parameters
n_replicas: 8
temperature:
  min: 300.0
  max: 400.0
temperature_method: "linear"
replex: 1000

# GROMACS settings
gromacs:
  gmx_mpi_command: "gmx_mpi"
  n_cpus: 4
  n_gpus: 1
  script_types: ["slurm", "localrun", "test"]

# Simulation settings
simulation:
  production_time: 100.0
  dt: 0.002
```

### Full Configuration Template (`configs/config_template.yaml`)
Contains all available options with detailed descriptions and examples.

## Usage

### Basic Usage
```bash
# Complete REST2 setup
python main.py -c configs/config_simple.yaml

# Custom output directory
python main.py -c configs/config_simple.yaml -o ./my_rest2_sim
```

### Advanced Usage
```bash
# Only validate configuration
python main.py -c configs/config_simple.yaml --validate-only

# Only generate execution scripts
python main.py -c configs/config_simple.yaml --scripts-only

# Verbose output
python main.py -c configs/config_simple.yaml --verbose
```

## Workflow Steps

The automation performs these steps:

1. **Structure Analysis**: Identifies target atoms and nearby residues
2. **Topology Merge**: Merges topology files using GROMACS tools
3. **Solute Selection**: Modifies topology file to mark REST2 atoms
4. **Replica Generation**: Creates replica directories and copies files
5. **Temperature Control**: Generates temperature-specific input files with PLUMED REST2
6. **Script Generation**: Creates execution scripts for your environment

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