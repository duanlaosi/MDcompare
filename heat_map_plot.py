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
                    raise ValueError(f"cannot read {filename}: cannot change '{parts[0]}' or '{parts[1]}' into float \n wrong line:{line.strip()}") from e

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


def create_groups(data, frame_number, group_number, data_name):
    """
    Split data into groups and print information about each group
    """
    grouped_data = []
    step = frame_number // group_number
    
    for i in range(group_number):
        start_idx = i * step
        # For the last group, include all remaining frames
        end_idx = start_idx + step if i < group_number - 1 else frame_number
        
        # Extract data into lists
        current_group = data[start_idx:end_idx]
        grouped_data.append(current_group)
        
        # Print current group information
        print(f"{data_name} Group {i+1}: from index {start_idx} to index {end_idx-1}")
        print(f"{data_name} group {i+1}: {len(current_group)} elements")
        
        # Print first and last element for verification
        if len(current_group) > 0:
            print(f"First element: {current_group[0]}, Last element: {current_group[-1]}")
        else:
            print("Warning: Empty group detected!")
    
    return grouped_data


def perform_t_tests(grouped_data1, grouped_data2, group_number, label1, label2):
    """
    Perform t-tests between groups of two datasets and return p-values
    """
    print(f"Performing t-tests between {label1} and {label2} groups...")
    
    # Define two-dimensional array to save p-values
    p_values = np.zeros((group_number, group_number))
    
    # T-test for every group and save p-value
    for i in range(group_number):
        for j in range(group_number):
            # Perform t-test between groups
            t_stat, p_value = stats.ttest_ind(grouped_data1[i], grouped_data2[j], equal_var=False)
            p_values[i, j] = p_value
            print(f"T-test between {label1} group {i+1} and {label2} group {j+1}: p-value = {p_value:.4f}")
    
    return p_values


def create_heatmap(p_values, group_number, xlabel, ylabel, title, save_path):
    """
    Create and save a heatmap of p-values
    """
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(p_values,
                    cmap='YlOrRd',
                    xticklabels=range(1, (group_number + 1)),
                    yticklabels=range(1, (group_number + 1)),
                    annot=True,
                    fmt=".2f")
    plt.xlabel(f'{xlabel} Group')
    plt.ylabel(f'{ylabel} Group')
    plt.title(title)

    # Configure colorbar
    cbar = ax.collections[0].colorbar
    cbar.set_label('p-value')

    # Save and display plot
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved heatmap to {save_path}")
    plt.close()


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

    # Create groups for PDB and AF data
    pdb_grouped_data = create_groups(pdb_projection_PC1, args.frame_number, args.group_number, "PDB")
    af_grouped_data = create_groups(af_projection_PC1, args.frame_number, args.group_number, "AF")
    
    # Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    # 1. Perform t-tests between AF and PDB groups (original functionality)
    p_values_af_vs_pdb = perform_t_tests(af_grouped_data, pdb_grouped_data, args.group_number, "AF", "PDB")
    
    # 2. Perform t-tests within PDB groups
    p_values_pdb_internal = perform_t_tests(pdb_grouped_data, pdb_grouped_data, args.group_number, "PDB", "PDB")
    
    # 3. Perform t-tests within AF groups
    p_values_af_internal = perform_t_tests(af_grouped_data, af_grouped_data, args.group_number, "AF", "AF")
    
    print("\nGenerating plots...")
    
    # Create and save heatmaps
    # 1. AF vs PDB (original)
    create_heatmap(
        p_values_af_vs_pdb, 
        args.group_number, 
        "PDB", 
        "AlphaFold", 
        f'P-values Between AF and PDB Groups', 
        os.path.join(args.save_dir, f'heatmap_{args.save_time}_{args.group_number}_af_vs_pdb.png')
    )
    
    # 2. PDB internal comparison
    create_heatmap(
        p_values_pdb_internal, 
        args.group_number, 
        "PDB", 
        "PDB", 
        f'P-values Within PDB Groups', 
        os.path.join(args.save_dir, f'heatmap_{args.save_time}_{args.group_number}_pdb_internal.png')
    )
    
    # 3. AF internal comparison
    create_heatmap(
        p_values_af_internal, 
        args.group_number, 
        "AlphaFold", 
        "AlphaFold", 
        f'P-values Within AF Groups', 
        os.path.join(args.save_dir, f'heatmap_{args.save_time}_{args.group_number}af_internal.png')
    )


if __name__ == "__main__":
    main()
