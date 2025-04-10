#!/usr/bin/env python3
"""
Script to update figure paths in Python files to reflect new folder structure.
"""

import os
import re
import glob

def update_file(file_path):
    """Update image paths in a Python file."""
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Update Monterey Bay figure paths
    updated_content = re.sub(
        r'(["\'])monterey_([^"\']*\.png)(["\'])', 
        r'\1figures/monterey/monterey_\2\3', 
        content
    )
    
    # Update general air-sea figure paths - wind, temperature, humidity effects
    updated_content = re.sub(
        r'(["\'])([^"\']*_effect\.png)(["\'])', 
        r'\1figures/general_airsea/\2\3', 
        updated_content
    )
    
    # Update general air-sea figure paths - comparisons
    updated_content = re.sub(
        r'(["\'])([^"\']*_comparison\.png)(["\'])', 
        r'\1figures/general_airsea/\2\3', 
        updated_content
    )
    
    # Update wind profile path
    updated_content = re.sub(
        r'(["\'])(wind_profile\.png)(["\'])', 
        r'\1figures/general_airsea/\2\3', 
        updated_content
    )
    
    # Check if any changes were made
    if content != updated_content:
        with open(file_path, 'w') as f:
            f.write(updated_content)
        return True
    return False

def main():
    """Update all Python files in the current directory."""
    python_files = glob.glob("*.py")
    
    updated_files = []
    for py_file in python_files:
        if py_file == "update_figure_paths.py":
            continue  # Skip this script
        
        if update_file(py_file):
            updated_files.append(py_file)
    
    print(f"Updated {len(updated_files)} files:")
    for file in updated_files:
        print(f" - {file}")

if __name__ == "__main__":
    main() 