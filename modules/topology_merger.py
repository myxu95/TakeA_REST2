#!/usr/bin/env python3
"""
REST2 Topology Merger Module
Simple topology merging using gmx grompp -pp
"""

import subprocess
from pathlib import Path


def merge_topology_files(main_topology: str, output_topology: str, 
                        structure_file: str) -> bool:
    """
    Merge topology files using gmx grompp -pp
    
    Args:
        main_topology: Main topology file (.top)
        output_topology: Output processed topology file
        structure_file: Structure file (.gro, .pdb)
        
    Returns:
        True if merge successful
    """
    try:
        # Check if processed.top already exists
        if Path(output_topology).exists():
            print(f"✓ Processed topology file already exists: {output_topology}")
            return True
        
        # Check input files
        main_top_path = Path(main_topology)
        structure_path = Path(structure_file)
        
        if not main_top_path.exists():
            print(f"✗ Main topology file not found: {main_topology}")
            return False
        
        if not structure_path.exists():
            print(f"✗ Structure file not found: {structure_file}")
            return False
        
        # Use template MDP file
        mdp_file = "templates/temp.mdp"
        mdp_path = Path(mdp_file)
        
        if not mdp_path.exists():
            print(f"✗ Template MDP file not found: {mdp_file}")
            return False
        
        print(f"Merging topology files...")
        print(f"  Main topology: {main_topology}")
        print(f"  Structure file: {structure_file}")
        print(f"  Template MDP: {mdp_file}")
        print(f"  Output: {output_topology}")
        
        # Prepare grompp command
        cmd = [
            'gmx', 'grompp',
            '-c', str(Path(structure_file).absolute()),
            '-p', str(Path(main_topology).absolute()),
            '-f', str(Path(mdp_file).absolute()),
            '-o', 'temp.tpr',
            '-pp', str(Path(output_topology).absolute())
        ]
        
        # Run grompp
        print(f"  Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=main_top_path.parent
        )
        
        if result.returncode != 0:
            print(f"  Error output: {result.stderr}")
            return False
        
        # Check if output file was created
        output_path = Path(output_topology)
        if output_path.exists():
            print(f"✓ Topology merged successfully: {output_topology}")
            
            # Clean up temporary files
            temp_tpr = Path('temp.tpr')
            if temp_tpr.exists():
                temp_tpr.unlink()
                print(f"✓ Cleaned up temporary file: temp.tpr")
            
            return True
        else:
            print("✗ Output topology file was not created")
            return False
            
    except Exception as e:
        print(f"✗ Topology merge failed: {e}")
        return False


def main():
    """Test topology merger"""
    # Test with example files
    main_top = "example/topol.top"
    structure = "example/md.gro"
    output_top = "example/processed.top"
    
    if Path(main_top).exists():
        print("Testing topology merger...")
        success = merge_topology_files(main_top, output_top, structure)
        
        if success:
            print("✓ Topology merger test passed")
        else:
            print("✗ Topology merger test failed")
    else:
        print("Example topology file not found, skipping test")


if __name__ == "__main__":
    main() 