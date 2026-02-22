#!/home/duanqi/miniconda3/envs/py36_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
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
    parser = argparse.ArgumentParser(description='Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_file', type=str,
                        help='Path to PDB XVG file')
    
    parser.add_argument('af_file', type=str,
                        help='Path to AF XVG file')
    
    parser.add_argument('frame_number', type=int,
                        help='Number of frames')
    
    parser.add_argument('group_number', type=int,
                        help='Number of groups')
    
    parser.add_argument('min_interval', type=float,
                        help='Minimum sample interval')
    
    parser.add_argument('max_interval', type=float,
                        help='Maximum sample interval')
    
    parser.add_argument('num_intervals', type=int,
                        help='Number of intervals')
    
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

def perform_tests(grouped_data1, grouped_data2, group_number, label1, label2):
    """
    Perform levene-tests and t-tests between groups of two datasets and return p-values
    """
    print(f"Performing levene-tests and t-tests between {label1} and {label2} groups...")
    
    # Define two-dimensional array to save p-values
    p_levenes = np.zeros((group_number, group_number)) 
    p_ttests = np.zeros((group_number, group_number))
    
    # levene-test and t-test for every group and save p-value
    for i in range(group_number):
        for j in range(group_number):
            # Perform t-test between groups
            levene_stat, p_levene = stats.levene(grouped_data1[i], grouped_data2[j])
            equal_var = True if p_levene > 0.05 else False
            t_stat, p_ttest = stats.ttest_ind(grouped_data1[i], grouped_data2[j], equal_var=equal_var)
            p_levenes[i, j] = p_levene
            p_ttests[i, j] = p_ttest
            print(f"Levene-test between {label1} group {i+1} and {label2} group {j+1}: p-levene = {p_levene:.4f}")
            print(f"T-test between {label1} group {i+1} and {label2} group {j+1}: p-ttest = {p_ttest:.4f}")
    
    return p_levenes, p_ttests

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
    args = parse_arguments()

    # Get data
    try:
        print("Reading PDB file...")
        pdb_time, pdb_RMSD = read_xvg(args.pdb_file)
        print("Reading AF file...")
        af_time, af_RMSD = read_xvg(args.af_file)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)

    #seperate interval using num_intervals to equally seperate
    min_int = int(args.min_interval * 10)
    max_int = int(args.max_interval * 10)
    int_intervals = np.linspace(min_int, max_int, args.num_intervals, dtype=int)
    sample_intervals = int_intervals / 10.0
    step_sizes = [max(1, int_val) for int_val in int_intervals]
    print("Sample intervals:", sample_intervals)
    print("Step sizes:", step_sizes)
    print("\nCalculating statistics...")
    # Loop through stepsize from 1 to 10(0.1ns sample interval to 1ns sample interval correspond to trajectory frame 1 to 10)
    for i, sample_interval in enumerate(sample_intervals):
        step = step_sizes[i]
        # Calculate the sampling for PDB data
        pdb_sample_RMSD = pdb_RMSD[::step]
        pdb_sample_time = pdb_time[::step]
        
        # Calculate the sampling for AF data
        af_sample_RMSD = af_RMSD[::step]
        af_sample_time = af_time[::step]

        # Calculate sampled frame number(equal to int(np.ceil(args.frame_number / step)), if pdb frame is equivalent to af frame )
        pdb_sample_frame_number = len(pdb_sample_RMSD)
        af_sample_frame_number = len(af_sample_RMSD)

        # Create groups for PDB and AF data
        pdb_grouped_data = create_groups(pdb_sample_RMSD, pdb_sample_frame_number, args.group_number, "PDB")
        af_grouped_data = create_groups(af_sample_RMSD, af_sample_frame_number, args.group_number, "AF")
    
        # 1. Perform t-tests between AF and PDB groups (original functionality)
        p_levenes_af_vs_pdb, p_ttests_af_vs_pdb = perform_tests(af_grouped_data, pdb_grouped_data, args.group_number, "AF", "PDB")
        
        # 2. Perform t-tests within PDB groups
        p_levenes_pdb_internal, p_ttests_pdb_internal = perform_tests(pdb_grouped_data, pdb_grouped_data, args.group_number, "PDB", "PDB")
        
        # 3. Perform t-tests within AF groups
        p_levenes_af_internal, p_ttests_af_internal = perform_tests(af_grouped_data, af_grouped_data, args.group_number, "AF", "AF")
        
        print("\nGenerating plots...")

        # Create and save heatmaps
        # 1. AF vs PDB (original)
        # 1.1 Levene test
        create_heatmap(
            p_levenes_af_vs_pdb, 
            args.group_number, 
            "PDB", 
            "AlphaFold", 
            f'Levene test P-values Between AF and PDB Groups', 
            os.path.join(args.save_dir, f'heatmap_levene_{args.save_time}_{args.group_number}_{step}_af_vs_pdb.png')
        )
        #1.2 T test
        create_heatmap(
            p_ttests_af_vs_pdb, 
            args.group_number, 
            "PDB", 
            "AlphaFold", 
            f'T test P-values Between AF and PDB Groups', 
            os.path.join(args.save_dir, f'heatmap_t_{args.save_time}_{args.group_number}_{step}_af_vs_pdb.png')
        )        

        # 2. PDB internal comparison
        # 2.1 Levene test
        create_heatmap(
            p_levenes_pdb_internal, 
            args.group_number, 
            "PDB", 
            "PDB", 
            f'Levene test P-values Within PDB Groups', 
            os.path.join(args.save_dir, f'heatmap_levene_{args.save_time}_{args.group_number}_{step}_pdb_internal.png')
        )
        # 2.2 T test
        create_heatmap(
            p_ttests_pdb_internal, 
            args.group_number, 
            "PDB", 
            "PDB", 
            f'T test P-values Within PDB Groups', 
            os.path.join(args.save_dir, f'heatmap_t_{args.save_time}_{args.group_number}_{step}_pdb_internal.png')
        )
        
        # 3. AF internal comparison
        # 3.1 Levene test
        create_heatmap(
            p_levenes_af_internal, 
            args.group_number, 
            "AlphaFold", 
            "AlphaFold", 
            f'Levene test P-values Within AF Groups', 
            os.path.join(args.save_dir, f'heatmap_levene_{args.save_time}_{args.group_number}_{step}_af_internal.png')
        )
        # 3.2 T test
        create_heatmap(
            p_ttests_af_internal, 
            args.group_number, 
            "AlphaFold", 
            "AlphaFold", 
            f'T test P-values Within AF Groups', 
            os.path.join(args.save_dir, f'heatmap_t_{args.save_time}_{args.group_number}_{step}_af_internal.png')
        )

if __name__ == "__main__":
    main()