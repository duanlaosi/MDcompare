#!/home/duanqi/miniconda3/envs/py36_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy import stats
import seaborn as sns
import os
import sys
import argparse

def read_xvg(filename):
    """
    read the data from xvg file and skip the annotation
    """
    x_values = []
    y_values = []
    
    with open(filename, 'r') as file:
        for line in file:

            # skip title and annotations
            if line.startswith('#') or line.startswith('@'):
                continue
            
            # split data lines
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    x_values.append(x)
                    y_values.append(y)
                except ValueError as e:
                    raise ValueError(f"cannot read traj_fit_eigenval.xvg:cannot change '{parts[0]}' or '{parts[1]}' into float \n wrong line:{line.strip()}") from e 
                    


        

    return np.array(x_values), np.array(y_values)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='PCA Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_file', type=str,
                        help='Path to PDB XVG file')
    
    parser.add_argument('af_file', type=str,
                        help='Path to AF XVG file')
    
    parser.add_argument('frame_number', type=int,
                        help='Number of frames')
    
    parser.add_argument('group_number', type=int,
                        help='Number of groups')
    
    parser.add_argument('save_dir', type=str,
                        help='Directory to save the plots')
    
    parser.add_argument('save_time', type=str,
                        help='Time to save the plots')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  PDB file: {args.pdb_file}")
    print(f"  AF file: {args.af_file}")
    print(f"  Number of frames: {args.frame_number}")
    print(f"  Number of groups: {args.group_number}")
    print(f"  Save directory: {args.save_dir}")
    
    # Check if the files exist
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file '{args.pdb_file}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(args.af_file):
        print(f"Error: AF file '{args.af_file}' does not exist.")
        sys.exit(1)


    # Get data
    try:
        print("Reading PDB file...")
        pdb_projection_PC1, pdb_projection_PC2 = read_xvg(args.pdb_file)
        print("Reading AF file...")
        af_projection_PC1, af_projection_PC2 = read_xvg(args.af_file)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)




    #define lists to seperate pdb_projection_PC1&af_projection_PC1 into groups
    pdb_grouped_data = []
    af_grouped_data = []
    # Calculate the step size based on the number of frames and groups
    step = args.frame_number // args.group_number
    # Group the data
    for i in range(args.group_number):
        start_idx = i * step
        # For the last group, include all remaining frames
        end_idx = start_idx + step if i < args.group_number - 1 else args.frame_number
        
        # Extract data into lists
        pdb_current_group = pdb_projection_PC1[start_idx:end_idx]
        af_current_group = af_projection_PC1[start_idx:end_idx]
        pdb_grouped_data.append(pdb_current_group)
        af_grouped_data.append(af_current_group)
        
        # Print current group information
        print(f"Group {i+1}: from index {start_idx} to index {end_idx-1}")
        print(f"PDB group {i+1}: {len(pdb_current_group)} elements")
        print(f"AF group {i+1}: {len(af_current_group)} elements")
        
        # Print first and last element for verification
        if len(pdb_current_group) > 0 and len(af_current_group) > 0:
            print(f"First element in PDB: {pdb_current_group[0]}, Last element: {pdb_current_group[-1]}")
            print(f"First element in AF: {af_current_group[0]}, Last element: {af_current_group[-1]}")
        else:
            print("Warning: Empty group detected!")
  
    # T-test & heat-map visualization
    print("Performing t-tests between groups...")

    # Define two-dimensional array to save p-values
    p_values = np.zeros((args.group_number, args.group_number))
    # T-test for every group and save p-value
    for i in range(args.group_number):
        for j in range(args.group_number):
            # Perform t-test between AF and PDB groups
            t_stat, p_value = stats.ttest_ind(af_grouped_data[i], pdb_grouped_data[j], equal_var=False)
            p_values[i, j] = p_value
            print(f"T-test between AF group {i+1} and PDB group {j+1}: p-value = {p_value:.4f}")


    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    print("\nGenerating plots...")
    # Create heatmap for -log10(p-values)
    ax = sns.heatmap(p_values,
                    cmap='YlOrRd',
                    xticklabels=range(1, (args.group_number + 1)),
                    yticklabels=range(1, (args.group_number + 1)),
                    annot=True,
                    fmt=".2f")
    plt.xlabel('PDB Group')
    plt.ylabel('AlphaFold Group')
    plt.title(f'P-value for Different Groups')

    #Configure colorbar
    cbar = ax.collections[0].colorbar
    cbar.set_label('p-value')

    # Save and display plot
    plt.tight_layout()
    save_path = os.path.join(args.save_dir, f'heatmap_{args.save_time}_t-pvalue.png')
    plt.savefig(save_path, dpi=300)
    print(f"Saved ttest heatmap to {save_path}")
if __name__ == "__main__":
    main()






