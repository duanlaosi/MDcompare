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
            #skip title and annotations
            if line.startswith('#') or line.startswith('@'):
                continue
            
            #split data lines
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
    
    #Define positional arguments (required, no default values)
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
    
    parser.add_argument('--start_interval', type=float, default=0.1,
                        help='Starting sample interval in ns (default: 0.1)')
    
    parser.add_argument('--end_interval', type=float, default=10.0,
                        help='Ending sample interval in ns (default: 10.0)')
    
    parser.add_argument('--interval_steps', type=int, default=10,
                        help='Number of interval steps between start and end (default: 10)')
    
    return parser.parse_args()


def create_groups(data, frame_number, group_number, data_name):
    """
    Split data into groups and print information about each group
    """
    grouped_data = []
    step = frame_number // group_number
    
    for i in range(group_number):
        start_idx = i * step
        #For the last group, include all remaining frames
        end_idx = start_idx + step if i < group_number - 1 else frame_number
        
        #Extract data into lists
        current_group = data[start_idx:end_idx]
        grouped_data.append(current_group)
        
        #Print current group information
        print(f"{data_name} Group {i+1}: from index {start_idx} to index {end_idx-1}")
        print(f"{data_name} group {i+1}: {len(current_group)} elements")
        
        #Print first and last element for verification
        if len(current_group) > 0:
            print(f"First element: {current_group[0]}, Last element: {current_group[-1]}")
        else:
            print("Warning: Empty group detected!")
    
    return grouped_data

def sample_data_with_interval(data, interval, original_interval=0.1):
    """
    Sample data with a specific interval (in ns)
    
    Parameters:
    - data: original data array
    - interval: desired sampling interval in ns
    - original_interval: original interval between frames in ns (default 0.1ns)
    
    Returns:
    - sampled data array
    """
    #Calculate the step size
    step = int(round(interval / original_interval))
    
    #Ensure step is at least 1
    step = max(1, step)
    
    #Sample the data
    sampled_data = data[::step]
    
    return sampled_data

def perform_statistical_tests(grouped_data1, grouped_data2, group_number, label1, label2):
    """
    Perform t-tests and Levene tests between groups of two datasets and return p-values
    """
    print(f"Performing statistical tests between {label1} and {label2} groups...")
    
    #Define two-dimensional arrays to save p-values
    t_test_p_values = np.zeros((group_number, group_number))
    levene_test_p_values = np.zeros((group_number, group_number))
    
    #Perform tests for every group and save p-values
    for i in range(group_number):
        for j in range(group_number):
            #Perform t-test between groups (for means comparison)
            t_stat, t_p_value = stats.ttest_ind(grouped_data1[i], grouped_data2[j], equal_var=False)
            t_test_p_values[i, j] = t_p_value
            
            #Perform Levene test between groups (for variance comparison)
            levene_stat, levene_p_value = stats.levene(grouped_data1[i], grouped_data2[j])
            levene_test_p_values[i, j] = levene_p_value
            
            print(f"Tests between {label1} group {i+1} and {label2} group {j+1}:")
            print(f"  T-test p-value = {t_p_value:.4f}")
            print(f"  Levene test p-value = {levene_p_value:.4f}")
    
    return t_test_p_values, levene_test_p_values


def create_heatmap(p_values, group_number, xlabel, ylabel, title, save_path, test_type="t-test"):
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
    plt.title(f'{title} - {test_type}')

    #Configure colorbar
    cbar = ax.collections[0].colorbar
    cbar.set_label('p-value')

    #Save and display plot
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Saved heatmap to {save_path}")
    plt.close()

def plot_p_values_vs_intervals(intervals, t_test_p_values_list, levene_test_p_values_list, 
                           comparison_type, save_dir, save_time):
    """
    Plot p-values across different sampling intervals
    
    Parameters:
    - intervals: list of sampling intervals
    - t_test_p_values_list: list of t-test p-values for each interval
    - levene_test_p_values_list: list of Levene test p-values for each interval
    - comparison_type: type of comparison (e.g., "AF_vs_PDB", "PDB_internal", "AF_internal")
    - save_dir: directory to save plots
    - save_time: timestamp for file naming
    """
    #Get the dimensions of p-values matrices
    n_rows, n_cols = t_test_p_values_list[0].shape
    
    #Create a figure with subplots for each group pair
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 15), sharex=True)

    #Plot for each group pair
    for i in range(n_rows):
        for j in range(n_cols):
            ax = axes[i, j] if n_rows > 1 else axes[j]
            
            #Extract p-values for this group pair across all intervals
            t_p_values = [p_vals[i, j] for p_vals in t_test_p_values_list]
            levene_p_values = [p_vals[i, j] for p_vals in levene_test_p_values_list]
            
            #Plot p-values
            ax.plot(intervals, t_p_values, 'bo-', label='T-test')
            ax.plot(intervals, levene_p_values, 'ro-', label='Levene test')
            
            #Add horizontal line at p=0.05 for significance threshold
            ax.axhline(y=0.05, color='k', linestyle='--', alpha=0.5)
            
            #Set title and labels
            if comparison_type == "AF_vs_PDB":
                ax.set_title(f'AF Group {i+1} vs PDB Group {j+1}')
            elif comparison_type == "PDB_internal":
                ax.set_title(f'PDB Group {i+1} vs PDB Group {j+1}')
            else:  #AF_internal
                ax.set_title(f'AF Group {i+1} vs AF Group {j+1}')
            
            #Only set x-label for bottom row
            if i == n_rows - 1:
                ax.set_xlabel('Sampling Interval (ns)')
            
            #Only set y-label for leftmost column
            if j == 0:
                ax.set_ylabel('p-value')
            
            #Set y-axis limits
            ax.set_ylim(0, 1.1)
            
            #Add legend only to the first subplot
            if i == 0 and j == 0:
                ax.legend()
    
    #Adjust layout
    plt.tight_layout()
    
    #Save the figure
    save_path = os.path.join(save_dir, f'interval_p_values_{save_time}_{comparison_type}.png')
    plt.savefig(save_path, dpi=300)
    print(f"Saved interval p-values plot to {save_path}")
    plt.close()


def main():
    #Parse command line arguments
    args = parse_arguments()
    
    #Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  PDB file: {args.pdb_file}")
    print(f"  AF file: {args.af_file}")
    print(f"  Number of frames: {args.frame_number}")
    print(f"  Number of groups: {args.group_number}")
    print(f"  Save directory: {args.save_dir}")
    print(f"  Sample interval range: {args.start_interval}ns to {args.end_interval}ns")
    print(f"  Number of interval steps: {args.interval_steps}")
    
    #Check if the files exist
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file '{args.pdb_file}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(args.af_file):
        print(f"Error: AF file '{args.af_file}' does not exist.")
        sys.exit(1)

    #Get data
    try:
        print("Reading PDB file...")
        pdb_projection_PC1, pdb_projection_PC2 = read_xvg(args.pdb_file)
        print("Reading AF file...")
        af_projection_PC1, af_projection_PC2 = read_xvg(args.af_file)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)
    
    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    #Generate the intervals to test
    intervals = np.linspace(args.start_interval, args.end_interval, args.interval_steps)
    print(f"Testing sample intervals: {intervals}")
    
    #Lists to store p-values for each interval
    t_test_p_values_af_vs_pdb_list = []
    levene_p_values_af_vs_pdb_list = []
    
    t_test_p_values_pdb_internal_list = []
    levene_p_values_pdb_internal_list = []
    
    t_test_p_values_af_internal_list = []
    levene_p_values_af_internal_list = []
    
    #For each sampling interval
    for interval in intervals:
        print(f"\n=== Processing sample interval: {interval}ns ===")
        
        #Sample the data for this interval
        pdb_sampled = sample_data_with_interval(pdb_projection_PC1, interval)
        af_sampled = sample_data_with_interval(af_projection_PC1, interval)
        
        #Calculate effective frame number after sampling
        effective_frame_number = min(len(pdb_sampled), len(af_sampled))
        print(f"Effective frame number after sampling: {effective_frame_number}")
        
        #Create groups for PDB and AF data
        pdb_grouped_data = create_groups(pdb_sampled[:effective_frame_number], 
                                         effective_frame_number, 
                                         args.group_number, 
                                         f"PDB (interval {interval}ns)")
        
        af_grouped_data = create_groups(af_sampled[:effective_frame_number], 
                                        effective_frame_number, 
                                        args.group_number, 
                                        f"AF (interval {interval}ns)")
        
        #1. Perform statistical tests between AF and PDB groups
        t_p_values_af_vs_pdb, levene_p_values_af_vs_pdb = perform_statistical_tests(
            af_grouped_data, pdb_grouped_data, args.group_number, "AF", "PDB")
        
        #2. Perform statistical tests within PDB groups
        t_p_values_pdb_internal, levene_p_values_pdb_internal = perform_statistical_tests(
            pdb_grouped_data, pdb_grouped_data, args.group_number, "PDB", "PDB")
        
        #3. Perform statistical tests within AF groups
        t_p_values_af_internal, levene_p_values_af_internal = perform_statistical_tests(
            af_grouped_data, af_grouped_data, args.group_number, "AF", "AF")
        
        #Store the p-values for this interval
        t_test_p_values_af_vs_pdb_list.append(t_p_values_af_vs_pdb)
        levene_p_values_af_vs_pdb_list.append(levene_p_values_af_vs_pdb)
        
        t_test_p_values_pdb_internal_list.append(t_p_values_pdb_internal)
        levene_p_values_pdb_internal_list.append(levene_p_values_pdb_internal)
        
        t_test_p_values_af_internal_list.append(t_p_values_af_internal)
        levene_p_values_af_internal_list.append(levene_p_values_af_internal)
        
        #Create and save heatmaps for this interval
        print("\nGenerating heatmaps for interval:", interval)
        
        #1. AF vs PDB
        #T-test heatmap
        create_heatmap(
            t_p_values_af_vs_pdb, 
            args.group_number, 
            "PDB", 
            "AlphaFold", 
            f'P-values Between AF and PDB Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_af_vs_pdb_ttest_{interval}ns.png'),
            "T-test"
        )
        
        #Levene test heatmap
        create_heatmap(
            levene_p_values_af_vs_pdb, 
            args.group_number, 
            "PDB", 
            "AlphaFold", 
            f'P-values Between AF and PDB Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_af_vs_pdb_levene_{interval}ns.png'),
            "Levene test"
        )
        
        #2. PDB internal
        #T-test heatmap
        create_heatmap(
            t_p_values_pdb_internal, 
            args.group_number, 
            "PDB", 
            "PDB", 
            f'P-values Within PDB Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_pdb_internal_ttest_{interval}ns.png'),
            "T-test"
        )
        
        #Levene test heatmap
        create_heatmap(
            levene_p_values_pdb_internal, 
            args.group_number, 
            "PDB", 
            "PDB", 
            f'P-values Within PDB Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_pdb_internal_levene_{interval}ns.png'),
            "Levene test"
        )
        
        #3. AF internal
        #T-test heatmap
        create_heatmap(
            t_p_values_af_internal, 
            args.group_number, 
            "AlphaFold", 
            "AlphaFold", 
            f'P-values Within AF Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_af_internal_ttest_{interval}ns.png'),
            "T-test"
        )
        
        #Levene test heatmap
        create_heatmap(
            levene_p_values_af_internal, 
            args.group_number, 
            "AlphaFold", 
            "AlphaFold", 
            f'P-values Within AF Groups (Interval {interval}ns)', 
            os.path.join(args.save_dir, f'heatmap_{args.save_time}_af_internal_levene_{interval}ns.png'),
            "Levene test"
        )
    
    #Plot p-values vs intervals for each group pair
    print("\nGenerating p-values vs intervals plots...")
    
    #1. AF vs PDB
    plot_p_values_vs_intervals(
        intervals,
        t_test_p_values_af_vs_pdb_list,
        levene_p_values_af_vs_pdb_list,
        "AF_vs_PDB",
        args.save_dir,
        args.save_time
    )
    
    #2. PDB internal
    plot_p_values_vs_intervals(
        intervals,
        t_test_p_values_pdb_internal_list,
        levene_p_values_pdb_internal_list,
        "PDB_internal",
        args.save_dir,
        args.save_time
    )
    
    #3. AF internal
    plot_p_values_vs_intervals(
        intervals,
        t_test_p_values_af_internal_list,
        levene_p_values_af_internal_list,
        "AF_internal",
        args.save_dir,
        args.save_time
    )


if __name__ == "__main__":
    main()