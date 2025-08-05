#!/usr/bin/env python3
"""
REST2 File Utilities Module
Unified file operations for REST2 enhanced sampling simulations
"""

import os
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
import yaml


class FileOperationError(Exception):
    """File operation error"""
    pass


class FileUtils:
    """
    Unified file utilities for REST2 simulations
    Handles all file operations across different modules
    """
    
    @staticmethod
    def ensure_directory(path: Union[str, Path], create_parents: bool = True) -> Path:
        """
        Ensure directory exists, create if necessary
        
        Args:
            path: Directory path
            create_parents: Whether to create parent directories
            
        Returns:
            Path object for the directory
        """
        path_obj = Path(path)
        if create_parents:
            path_obj.mkdir(parents=True, exist_ok=True)
        else:
            path_obj.mkdir(exist_ok=True)
        return path_obj
    
    @staticmethod
    def safe_copy(src: Union[str, Path], dst: Union[str, Path], 
                  overwrite: bool = False, backup: bool = True) -> Path:
        """
        Safely copy file with optional backup
        
        Args:
            src: Source file path
            dst: Destination file path
            overwrite: Whether to overwrite existing file
            backup: Whether to create backup of existing file
            
        Returns:
            Path object for the destination file
            
        Raises:
            FileOperationError: If copy operation fails
        """
        src_path = Path(src)
        dst_path = Path(dst)
        
        if not src_path.exists():
            raise FileOperationError(f"Source file not found: {src_path}")
        
        # Create destination directory if needed
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Handle existing destination file
        if dst_path.exists():
            if not overwrite:
                raise FileOperationError(f"Destination file exists: {dst_path}")
            
            if backup:
                backup_path = dst_path.with_suffix(dst_path.suffix + '.backup')
                if not backup_path.exists():
                    shutil.copy2(dst_path, backup_path)
        
        # Copy file
        try:
            shutil.copy2(src_path, dst_path)
            return dst_path
        except Exception as e:
            raise FileOperationError(f"Failed to copy {src_path} to {dst_path}: {e}")
    
    @staticmethod
    def safe_write(content: str, file_path: Union[str, Path], 
                   overwrite: bool = False, backup: bool = True) -> Path:
        """
        Safely write content to file
        
        Args:
            content: Content to write
            file_path: Target file path
            overwrite: Whether to overwrite existing file
            backup: Whether to create backup of existing file
            
        Returns:
            Path object for the written file
        """
        path_obj = Path(file_path)
        
        # Create directory if needed
        path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        # Handle existing file
        if path_obj.exists():
            if not overwrite:
                raise FileOperationError(f"File exists: {path_obj}")
            
            if backup:
                backup_path = path_obj.with_suffix(path_obj.suffix + '.backup')
                if not backup_path.exists():
                    shutil.copy2(path_obj, backup_path)
        
        # Write content
        try:
            with open(path_obj, 'w', encoding='utf-8') as f:
                f.write(content)
            return path_obj
        except Exception as e:
            raise FileOperationError(f"Failed to write to {path_obj}: {e}")
    
    @staticmethod
    def safe_read(file_path: Union[str, Path]) -> str:
        """
        Safely read file content
        
        Args:
            file_path: File path to read
            
        Returns:
            File content as string
            
        Raises:
            FileOperationError: If read operation fails
        """
        path_obj = Path(file_path)
        
        if not path_obj.exists():
            raise FileOperationError(f"File not found: {path_obj}")
        
        try:
            with open(path_obj, 'r', encoding='utf-8') as f:
                return f.read()
        except Exception as e:
            raise FileOperationError(f"Failed to read {path_obj}: {e}")
    
    @staticmethod
    def load_yaml(file_path: Union[str, Path]) -> Dict[str, Any]:
        """
        Load YAML file
        
        Args:
            file_path: YAML file path
            
        Returns:
            Dictionary with YAML content
            
        Raises:
            FileOperationError: If YAML loading fails
        """
        try:
            content = FileUtils.safe_read(file_path)
            return yaml.safe_load(content) or {}
        except yaml.YAMLError as e:
            raise FileOperationError(f"Invalid YAML format in {file_path}: {e}")
        except Exception as e:
            raise FileOperationError(f"Failed to load YAML from {file_path}: {e}")
    
    @staticmethod
    def save_yaml(data: Dict[str, Any], file_path: Union[str, Path], 
                  default_flow_style: bool = False, indent: int = 2) -> Path:
        """
        Save data to YAML file
        
        Args:
            data: Data to save
            file_path: Target YAML file path
            default_flow_style: YAML flow style
            indent: YAML indentation
            
        Returns:
            Path object for the saved file
        """
        try:
            content = yaml.dump(data, default_flow_style=default_flow_style, indent=indent)
            return FileUtils.safe_write(content, file_path)
        except Exception as e:
            raise FileOperationError(f"Failed to save YAML to {file_path}: {e}")
    
    @staticmethod
    def find_files(directory: Union[str, Path], patterns: List[str], 
                   recursive: bool = True) -> List[Path]:
        """
        Find files matching patterns in directory
        
        Args:
            directory: Directory to search
            patterns: List of file patterns (glob)
            recursive: Whether to search recursively
            
        Returns:
            List of matching file paths
        """
        dir_path = Path(directory)
        if not dir_path.exists():
            return []
        
        found_files = []
        for pattern in patterns:
            if recursive:
                found_files.extend(dir_path.rglob(pattern))
            else:
                found_files.extend(dir_path.glob(pattern))
        
        return list(set(found_files))  # Remove duplicates
    
    @staticmethod
    def auto_detect_files(directory: Union[str, Path], 
                         file_mapping: Dict[str, List[str]]) -> Dict[str, Optional[str]]:
        """
        Auto-detect files in directory based on common names
        
        Args:
            directory: Directory to search
            file_mapping: Mapping of config keys to possible filenames
            
        Returns:
            Dictionary mapping config keys to found file paths
        """
        dir_path = Path(directory)
        if not dir_path.exists():
            return {key: None for key in file_mapping.keys()}
        
        results = {}
        for config_key, possible_names in file_mapping.items():
            results[config_key] = None
            for filename in possible_names:
                file_path = dir_path / filename
                if file_path.exists():
                    results[config_key] = str(file_path)
                    break
        
        return results
    
    @staticmethod
    def create_directory_structure(base_dir: Union[str, Path], 
                                 structure: Dict[str, Any]) -> Dict[str, Path]:
        """
        Create directory structure
        
        Args:
            base_dir: Base directory
            structure: Directory structure definition
            
        Returns:
            Dictionary mapping structure keys to created paths
        """
        base_path = Path(base_dir)
        created_paths = {}
        
        def create_structure_recursive(parent_path: Path, struct: Dict[str, Any], 
                                     path_dict: Dict[str, Path]) -> None:
            for key, value in struct.items():
                if isinstance(value, dict):
                    # Directory
                    dir_path = parent_path / key
                    dir_path.mkdir(parents=True, exist_ok=True)
                    path_dict[key] = dir_path
                    create_structure_recursive(dir_path, value, path_dict)
                else:
                    # File or other value
                    path_dict[key] = parent_path / value
        
        create_structure_recursive(base_path, structure, created_paths)
        return created_paths
    
    @staticmethod
    def validate_file_structure(directory: Union[str, Path], 
                              required_files: List[str], 
                              optional_files: List[str] = None) -> Dict[str, bool]:
        """
        Validate file structure in directory
        
        Args:
            directory: Directory to validate
            required_files: List of required files
            optional_files: List of optional files
            
        Returns:
            Dictionary mapping file names to existence status
        """
        dir_path = Path(directory)
        if not dir_path.exists():
            return {file: False for file in required_files + (optional_files or [])}
        
        results = {}
        
        # Check required files
        for filename in required_files:
            file_path = dir_path / filename
            results[filename] = file_path.exists()
        
        # Check optional files
        if optional_files:
            for filename in optional_files:
                file_path = dir_path / filename
                results[filename] = file_path.exists()
        
        return results
    
    @staticmethod
    def get_file_info(file_path: Union[str, Path]) -> Dict[str, Any]:
        """
        Get detailed file information
        
        Args:
            file_path: File path
            
        Returns:
            Dictionary with file information
        """
        path_obj = Path(file_path)
        
        if not path_obj.exists():
            return {
                'exists': False,
                'path': str(path_obj),
                'error': 'File not found'
            }
        
        try:
            stat = path_obj.stat()
            return {
                'exists': True,
                'path': str(path_obj),
                'size': stat.st_size,
                'modified': stat.st_mtime,
                'is_file': path_obj.is_file(),
                'is_dir': path_obj.is_dir(),
                'extension': path_obj.suffix,
                'name': path_obj.name,
                'parent': str(path_obj.parent)
            }
        except Exception as e:
            return {
                'exists': False,
                'path': str(path_obj),
                'error': str(e)
            }


def main():
    """Test file utilities"""
    try:
        print("Testing File Utils")
        print("=" * 40)
        
        # Test directory creation
        test_dir = Path("./test_file_utils")
        FileUtils.ensure_directory(test_dir)
        print(f"✓ Created directory: {test_dir}")
        
        # Test file writing
        test_file = test_dir / "test.txt"
        content = "Hello, File Utils!"
        FileUtils.safe_write(content, test_file)
        print(f"✓ Wrote file: {test_file}")
        
        # Test file reading
        read_content = FileUtils.safe_read(test_file)
        print(f"✓ Read content: {read_content}")
        
        # Test YAML operations
        yaml_data = {
            'test': 'value',
            'numbers': [1, 2, 3],
            'nested': {'key': 'value'}
        }
        yaml_file = test_dir / "test.yaml"
        FileUtils.save_yaml(yaml_data, yaml_file)
        print(f"✓ Saved YAML: {yaml_file}")
        
        loaded_data = FileUtils.load_yaml(yaml_file)
        print(f"✓ Loaded YAML: {loaded_data}")
        
        # Test file detection
        file_mapping = {
            'config': ['config.yaml', 'config.yml'],
            'data': ['data.txt', 'input.txt']
        }
        detected = FileUtils.auto_detect_files(test_dir, file_mapping)
        print(f"✓ Detected files: {detected}")
        
        # Test file validation
        validation = FileUtils.validate_file_structure(
            test_dir, 
            ['test.txt'], 
            ['test.yaml']
        )
        print(f"✓ File validation: {validation}")
        
        # Test file info
        file_info = FileUtils.get_file_info(test_file)
        print(f"✓ File info: {file_info}")
        
        # Cleanup
        shutil.rmtree(test_dir)
        print("✓ Cleanup completed")
        
        print("\n✓ File utils tests passed!")
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 