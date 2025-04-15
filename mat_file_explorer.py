import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

def explore_mat_file(file_path):
    """
    Explore the contents of a .mat file and print its structure
    """
    try:
        # Load the .mat file
        mat_contents = scipy.io.loadmat(file_path)
        
        print(f"\nExploring: {file_path}")
        print("-" * 50)
        
        # Print all variables in the file
        print("\nVariables found in the file:")
        for key in mat_contents.keys():
            if key not in ['__header__', '__version__', '__globals__']:
                data = mat_contents[key]
                print(f"\nVariable: {key}")
                print(f"Type: {type(data)}")
                print(f"Shape: {data.shape}")
                print(f"Data type: {data.dtype}")
                
                # If it's a small array, print its contents
                if data.size < 10:
                    print(f"Contents: {data}")
                
        return mat_contents
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def visualize_data(mat_contents, variable_name):
    """
    Create basic visualizations for the specified variable
    """
    if variable_name not in mat_contents:
        print(f"Variable {variable_name} not found in the file")
        return
    
    data = mat_contents[variable_name]
    
    # Create a figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle(f'Visualizations for {variable_name}')
    
    # 1D data visualization
    if len(data.shape) == 1:
        axes[0, 0].plot(data)
        axes[0, 0].set_title('Line Plot')
        axes[0, 0].set_xlabel('Index')
        axes[0, 0].set_ylabel('Value')
        
        axes[0, 1].hist(data.flatten(), bins=30)
        axes[0, 1].set_title('Histogram')
        axes[0, 1].set_xlabel('Value')
        axes[0, 1].set_ylabel('Count')
        
        # Hide unused subplots
        axes[1, 0].set_visible(False)
        axes[1, 1].set_visible(False)
    
    # 2D data visualization
    elif len(data.shape) == 2:
        im = axes[0, 0].imshow(data, cmap='viridis')
        axes[0, 0].set_title('2D Data')
        plt.colorbar(im, ax=axes[0, 0])
        
        axes[0, 1].hist(data.flatten(), bins=30)
        axes[0, 1].set_title('Histogram')
        axes[0, 1].set_xlabel('Value')
        axes[0, 1].set_ylabel('Count')
        
        # Plot mean along each dimension
        axes[1, 0].plot(np.mean(data, axis=0))
        axes[1, 0].set_title('Mean along first dimension')
        axes[1, 0].set_xlabel('Index')
        axes[1, 0].set_ylabel('Mean Value')
        
        axes[1, 1].plot(np.mean(data, axis=1))
        axes[1, 1].set_title('Mean along second dimension')
        axes[1, 1].set_xlabel('Index')
        axes[1, 1].set_ylabel('Mean Value')
    
    plt.tight_layout()
    plt.show()

def main():
    # Get the .mat file path from user
    file_path = input("Enter the path to your .mat file: ")
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    # Explore the file
    mat_contents = explore_mat_file(file_path)
    
    if mat_contents is not None:
        # Ask user which variable to visualize
        print("\nAvailable variables for visualization:")
        for key in mat_contents.keys():
            if key not in ['__header__', '__version__', '__globals__']:
                print(f"- {key}")
        
        variable_name = input("\nEnter the variable name to visualize (or press Enter to skip): ")
        if variable_name:
            visualize_data(mat_contents, variable_name)

if __name__ == "__main__":
    main() 