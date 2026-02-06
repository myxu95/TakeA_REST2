# REST2 Enhanced Sampling Automation

**Automated setup and execution of REST2 (Replica Exchange with Solute Tempering) enhanced sampling simulations using GROMACS and PLUMED.**

REST2 is an advanced enhanced sampling technique that selectively heats the solute (region of interest) while keeping the solvent at a reference temperature, providing more efficient conformational sampling compared to traditional replica exchange methods.

---

## Features

- ✅ **Automated REST2 Setup**: Complete workflow automation from structure analysis to script generation
- ✅ **Exponential Temperature Scaling**: Optimized replica distribution with more sampling at lower temperatures
- ✅ **Flexible Target Selection**: Automatic selection of solute and nearby residues using MDAnalysis
- ✅ **PLUMED Integration**: Automatic generation of PLUMED PARTIAL_TEMPERING configurations
- ✅ **Script Generation**: Create SLURM, local run, and test scripts automatically
- ✅ **Validation Framework**: Built-in validation for all configuration parameters
- ✅ **Trajectory-based Selection**: Optional contact analysis using MD trajectories

---

## Requirements

### Software Dependencies

- **Python**: 3.7+ (tested with 3.12)
- **GROMACS**: 2019+ with MPI support (tested with 2023.2)
- **PLUMED**: 2.5+ compiled with GROMACS (tested with 2.9.0)

### Python Packages

Install via conda (recommended):
```bash
conda env create -f environment.yml
conda activate rest2_env
```

Or install via pip:
```bash
pip install -r requirements.txt
```

Required packages:
- `numpy>=1.21.0`
- `scipy>=1.7.0`
- `MDAnalysis>=2.0.0`
- `PyYAML>=6.0`

---

## Installation

### Method 1: Using Conda (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd REST2_Project

# Create conda environment
conda env create -f environment.yml
conda activate rest2_env

# Verify installation
python main.py --help
```

### Method 2: Using Pip

```bash
# Clone the repository
git clone <repository-url>
cd REST2_Project

# Install dependencies
pip install -r requirements.txt

# Verify installation
python main.py --help
```

### Verify GROMACS and PLUMED

```bash
# Check GROMACS
gmx --version

# Check PLUMED
plumed --version

# Check if PLUMED is patched with GROMACS
gmx mdrun -h 2>&1 | grep -i plumed
```

---

## Quick Start

### 1. Prepare Your MD Results

Complete a standard MD simulation workflow:
```
Energy Minimization (EM) → NVT → NPT → Production MD
```

Ensure you have these files:
- `md.tpr` - Production MD run file
- `topol.top` - System topology file
- `md.gro` - Final structure file
- `md.xtc` (optional) - Trajectory for contact analysis

### 2. Configure Your System

Edit a configuration file (or use the example):

```yaml
# configs/my_system.yaml
target_type: peptide
T_min: 300.0              # Reference temperature (K)
T_max: 340.0              # Maximum temperature (K)
n_replicas: 8             # Number of replicas
scaling_method: exponential  # Exponential temperature distribution

# File paths
input_tpr: path/to/md.tpr
topology: path/to/topol.top
md_results_dir: path/to/MD_results
output_dir: ./rest2_simulation

# Target selection (MDAnalysis syntax)
target_selection: "protein"  # or "chainid A" or "resname LIG"
distance_range: 5.0          # Angstroms

# GROMACS settings
gromacs:
  gmx_mpi_command: gmx       # or gmx_mpi if MPI version available
  n_cpus: 8
  n_gpus: 2
```

### 3. Run the Setup

```bash
# Activate environment
conda activate rest2_env

# Validate configuration
python main.py -c configs/my_system.yaml --validate-only

# Run complete setup
python main.py -c configs/my_system.yaml

# Or with verbose output
python main.py -c configs/my_system.yaml --verbose
```

### 4. Execute REST2 Simulation

```bash
cd rest2_simulation

# Test the setup (short run)
./test_rest2.sh

# Submit to cluster
sbatch run_rest2.slurm

# Or run locally
./run_rest2_local.sh
```

---

## Temperature Scaling Method

This tool uses **exponential temperature scaling** for optimal replica distribution:

### Formula

```
T[i] = T_min × exp(i × log(T_max/T_min) / (n-1))
```

Or equivalently:
```
ratio = (T_max/T_min)^(1/(n-1))
T[i] = T_min × ratio^i
```

### Characteristics

- ✅ **More replicas at lower temperatures**: ΔT increases from low to high temperature
- ✅ **Better sampling efficiency**: Enhanced sampling where it matters most
- ✅ **Optimal for REST2**: Better exchange acceptance rates at low temperatures

### Example (300-340 K, 8 replicas)

```
Replica  Temperature  λ         ΔT
0        300.0 K      1.000000  5.41 K
1        305.4 K      0.982278  5.51 K
2        310.9 K      0.964871  5.61 K
3        316.5 K      0.947772  5.71 K
4        322.2 K      0.930976  5.81 K
5        328.1 K      0.914478  5.92 K
6        333.9 K      0.898272  6.03 K
7        340.0 K      0.882353  -
```

**Key advantage**: ΔT increases from 5.41K to 6.03K, providing denser replica coverage at lower temperatures where conformational barriers are typically higher.

---

## Configuration Options

### Basic Settings

```yaml
# Target specification
target_type: peptide          # peptide or small_molecule
target_selection: "chainid A" # MDAnalysis selection syntax

# Temperature settings
T_min: 300.0                  # Reference temperature (K)
T_max: 340.0                  # Maximum temperature (K)
n_replicas: 8                 # Number of replicas
scaling_method: exponential   # Exponential temperature distribution
replex: 1000                  # Exchange attempt interval (steps)

# Solute selection
distance_range: 5.0           # Cutoff distance (Angstroms)
use_trajectory: false         # Use MD trajectory for contact analysis
occupancy_threshold: 0.5      # Contact occupancy threshold (if using trajectory)
```

### File Settings

```yaml
# Input files
input_tpr: example/md.tpr
topology: example/topol.top
md_results_dir: example
plumed_dat: templates/plumed.dat  # Optional custom PLUMED template

# Output settings
output_dir: ./rest2_simulation
force_overwrite: true
```

### GROMACS Settings

```yaml
gromacs:
  gmx_mpi_command: gmx         # GROMACS command (gmx or gmx_mpi)
  n_cpus: 8                    # Total CPU cores
  n_gpus: 2                    # Total GPUs
  script_types:                # Scripts to generate
    - slurm                    # SLURM cluster submission
    - localrun                 # Local execution
    - test                     # Quick validation
```

### Target Selection Examples

```yaml
# Select specific chain
target_selection: "chainid A"

# Select ligand
target_selection: "resname LIG"

# Select residue range
target_selection: "resid 1-100"

# Select protein
target_selection: "protein"

# Complex selection
target_selection: "protein and (resid 1-50 or resid 100-150)"
```

For complete options, see `configs/config_template.yaml`

---

## Troubleshooting

### Common Issues

#### 1. Configuration Validation Failed

**Problem**: Configuration file has invalid parameters

**Solution**:
```bash
# Check configuration syntax
python main.py -c config.yaml --validate-only

# Common issues:
# - Missing required files
# - T_max <= T_min
# - Invalid target_selection syntax
# - Wrong file paths
```

#### 2. ModuleNotFoundError: No module named 'MDAnalysis'

**Problem**: Python dependencies not installed

**Solution**:
```bash
# If using conda
conda activate rest2_env
conda install -c conda-forge mdanalysis scipy

# If using pip
pip install MDAnalysis scipy numpy PyYAML
```

#### 3. GROMACS Command Not Found

**Problem**: `gmx` or `gmx_mpi` not in PATH

**Solution**:
```bash
# Check GROMACS installation
which gmx
which gmx_mpi

# If not found, load module (on cluster)
module load gromacs

# Update config file
gromacs:
  gmx_mpi_command: gmx  # or gmx_mpi
```

#### 4. No Atoms Found for Selection

**Problem**: Target selection doesn't match any atoms

**Solution**:
```bash
# Test selection manually:
python -c "import MDAnalysis as mda; u = mda.Universe('example/md.tpr'); print(u.select_atoms('chainid C'))"

# Try simpler selection
target_selection: "protein"  # Instead of specific chain
```

### Getting Help

For detailed error messages, use verbose mode:
```bash
python main.py -c config.yaml --verbose
```

---

## Running the Example

The project includes a complete working example:

```bash
# 1. Check example files
ls example/
# Output: md.gro, md.tpr, topol.top, charmm36-jul2021.ff/

# 2. Validate example configuration
python main.py -c configs/example_config.yaml --validate-only

# 3. Run the example setup
python main.py -c configs/example_config.yaml

# 4. Verify generated files
ls rest2_simulation/
ls rest2_simulation/replica_0/input/

# 5. Run quick test
cd rest2_simulation
./test_rest2.sh
```

---

## Generated Output Structure

After running the setup, you will get:

```
rest2_simulation/
├── replica_0/
│   ├── input/
│   │   ├── input.tpr               # GROMACS run file
│   │   ├── topol.top               # Scaled topology
│   │   ├── rest2.mdp               # MD parameters
│   │   ├── plumed.dat              # PLUMED configuration
│   │   └── index.ndx               # Index file
│   ├── output/                     # Simulation outputs (created during run)
│   └── replica_info.txt            # Replica information
├── replica_1/
├── ... (replica_2 to replica_7)
│
├── processed.top                   # Merged topology
├── rest2_topol.top                # REST2-modified topology
├── selected_residues.txt          # Selected residues info
├── temperature_summary.txt        # Temperature ladder info
│
├── run_rest2.slurm                # SLURM submission script
├── run_rest2_local.sh             # Local execution script
└── test_rest2.sh                  # Quick test script
```

---

## Performance Tips

1. **GPU Allocation**: Use 1 GPU per replica for optimal performance
   ```yaml
   n_gpus: 8  # For 8 replicas
   ```

2. **CPU Allocation**: Match CPUs to replicas
   ```yaml
   n_cpus: 8  # For 8 replicas
   ```

3. **Exchange Frequency**: Start with conservative values
   ```yaml
   replex: 1000  # Exchange every 1000 steps (2 ps with dt=0.002)
   ```

4. **Number of Replicas**: More replicas = better exchange, but higher cost
   - Small systems: 4-8 replicas
   - Large systems: 8-16 replicas
   - Rule of thumb: Exchange acceptance ~20-30%

5. **Temperature Range**: Keep T_max/T_min < 1.3 for good exchange
   ```yaml
   T_min: 300.0
   T_max: 380.0  # Ratio = 1.27 ✓
   ```

---

## Command-Line Options

```bash
# Main options
python main.py -c CONFIG [-o OUTPUT_DIR] [OPTIONS]

Options:
  -c, --config CONFIG       Configuration file (required)
  -o, --output-dir DIR      Override output directory
  --validate-only           Only validate configuration
  --scripts-only            Only generate execution scripts
  --verbose, -v             Verbose output
  -h, --help                Show help message
```

---

## References

### REST2 Method
- Wang, L., Friesner, R. A., & Berne, B. J. (2011). *Replica exchange with solute scaling: a more efficient version of replica exchange with solute tempering (REST2)*. The Journal of Physical Chemistry B, 115(30), 9431-9438.

### PLUMED
- Tribello, G. A., Bonomi, M., Branduardi, D., Camilloni, C., & Bussi, G. (2014). *PLUMED 2: New feathers for an old bird*. Computer Physics Communications, 185(2), 604-613.
- PLUMED REST2 Documentation: https://www.plumed.org/doc-v2.9/user-doc/html/_p_a_r_t_i_a_l__t_e_m_p_e_r_i_n_g.html

### GROMACS
- Abraham, M. J., et al. (2015). *GROMACS: High performance molecular simulations*. SoftwareX, 1, 19-25.
- GROMACS Manual: https://manual.gromacs.org/

### MDAnalysis
- Michaud-Agrawal, N., et al. (2011). *MDAnalysis: a toolkit for the analysis of molecular dynamics simulations*. Journal of computational chemistry, 32(10), 2319-2327.

---

## Citation

If you use this tool in your research, please cite the REST2 original paper:

```bibtex
@article{wang2011rest2,
  title={Replica exchange with solute scaling: a more efficient version of replica exchange with solute tempering (REST2)},
  author={Wang, Lingle and Friesner, Richard A and Berne, Bruce J},
  journal={The Journal of Physical Chemistry B},
  volume={115},
  number={30},
  pages={9431--9438},
  year={2011},
  publisher={ACS Publications}
}
```

---

## License

See [LICENSE](LICENSE) file for details.

---

## Changelog

### Version 2.0 (Current)
- ✨ Optimized for exponential temperature scaling
- ✨ Enhanced validation framework
- ✨ Better error messages and logging
- 📝 Comprehensive documentation
- 🐛 Bug fixes and performance improvements

### Version 1.0
- 🎉 Initial release
- ✅ Basic REST2 setup automation
- ✅ PLUMED integration
- ✅ Script generation

---

**Happy Sampling! 🚀**
