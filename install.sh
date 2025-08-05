#!/bin/bash
# REST2 Enhanced Sampling Project Installation Script

set -e

echo "=========================================="
echo "REST2 Enhanced Sampling Project Setup"
echo "=========================================="

# Check if Python 3.8+ is available
python_version=$(python3 --version 2>&1 | awk '{print $2}' | cut -d. -f1,2)
required_version="3.8"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then
    echo "Error: Python 3.8 or higher is required. Found: $python_version"
    exit 1
fi

echo "Python version: $python_version"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "Error: pip3 is not installed"
    exit 1
fi

echo "pip3 is available"

# Create virtual environment if not exists
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install Python dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt

# Check if GROMACS is available
if command -v gmx &> /dev/null; then
    gmx_version=$(gmx --version 2>&1 | head -n1)
    echo "GROMACS found: $gmx_version"
else
    echo "Warning: GROMACS not found in PATH"
    echo "   Please install GROMACS or ensure it's in your PATH"
fi

# Check if PLUMED is available
if command -v plumed &> /dev/null; then
    plumed_version=$(plumed --version 2>&1 | head -n1)
    echo "PLUMED found: $plumed_version"
else
    echo "Warning: PLUMED not found in PATH"
    echo "   Please install PLUMED or ensure it's in your PATH"
fi

# Make main script executable
chmod +x main.py

echo ""
echo "=========================================="
echo "Installation completed successfully!"
echo "=========================================="
echo ""
echo "To activate the environment:"
echo "  source venv/bin/activate"
echo ""
echo "To run the REST2 setup:"
echo "  python main.py -c configs/config_simple.yaml"
echo ""
echo "For help:"
echo "  python main.py --help"
echo "" 