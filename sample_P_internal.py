#!/home/duanqi/miniconda3/envs/py38_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
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

def sample_data_with_interval(pdb_data, af_data, min_interval, max_interval, num_intervals):
    """
    seperate interval using num_intervals to equally seperate

    Returns:
        sampled_data: Nested dictionaries in the format:
            {
                'PDB': {
                    interval: {'data': ..., 'frame_count': ...}
                },
                'AF': {
                    interval: {'data': ..., 'frame_count': ...}
                }
            }
    """
    pdb_sampled_data = {}
    af_sampled_data = {}
    sampled_data = {}

    print("Sampling data with intervals...")
    #Calculate the step size
    min_int = int(min_interval * 10)
    max_int = int(max_interval * 10)
    step = (max_int - min_int) // (num_intervals - 1)
    int_intervals = list(range(min_int, max_int + 1, step))
    sample_intervals = [val / 10 for val in int_intervals]
    step_sizes = [max(1, int_val) for int_val in int_intervals]
    print("Sample intervals:", sample_intervals)
    print("Step sizes:", step_sizes)
    print("\nCalculating statistics...")
    
    #Sample the data
    # Loop through stepsize from 1 to 10(0.1ns sample interval to 1ns sample interval correspond to trajectory frame 1 to 10)
    for i, sample_interval in enumerate(sample_intervals):
        step = step_sizes[i]
        # Calculate the sampling for PDB data
        pdb_sample = pdb_data[::step]
        pdb_sampled_data[sample_interval] = {
            'data': pdb_sample,
            'frame_count': len(pdb_sample)
        }
        
        # Calculate the sampling for AF data
        af_sample = af_data[::step]
        af_sampled_data[sample_interval] = {
            'data': af_sample,
            'frame_count': len(af_sample)
        }
    
        print(f"sample interval {sample_interval} ns: step={step}, PDB frame Number={len(pdb_sample)}, AF frame Number={len(af_sample)}")
    
    sampled_data= {
        'PDB': pdb_sampled_data,
        'AF': af_sampled_data
    }
    return sampled_data

def group_sampled_data(sampled_data, group_number, data_name):
    """
    Grouping data with different sampling intervals
    
    Parameters:
    sampled_data: Sampling data dictionary
    group_number: Number of groups
    data_name: Data type identifier (such as "PDB" or "AF")
    
    Returns:
    grouped_data: Nested dictionary in the format of {sample_interval: {group_id: group_data}}
    """
    grouped_data = {}
    
    # Use get method to safely access the data_name('PDB' or 'AF') key's data
    data_by_interval = sampled_data.get(data_name, {})

    for sample_interval, data in data_by_interval.items():
        frame_count = data['frame_count']
        frames_per_group = frame_count // group_number
        
        grouped_data[sample_interval] = {}
     
        for group_id in range(group_number):
            start_idx = group_id * frames_per_group
            # For the last group, include all remaining frames
            end_idx = start_idx + frames_per_group if group_id < group_number - 1 else frame_count
            
            # Extract data into lists
            group_data = data['data'][start_idx:end_idx]
            grouped_data[sample_interval][group_id] = group_data
            
            
            # Print group data for verification
            if len(group_data) > 0:
                print(f"{data_name} interval {sample_interval}ns group {group_id+1}: {len(group_data)}frame "
                f"(index {start_idx}-{end_idx-1})")
            else:
                print("Warning: Empty group detected!")
                raise ValueError(f"Empty group detected for {data_name} at interval {sample_interval}ns, group {group_id+1}. ")

    return grouped_data



def perform_tests_internal(grouped_data, group_number, label):
    """
    Perform Levene's test and t-test on grouped data at multiple sampling intervals, returning a dictionary of p-values

    Parameters:
        grouped_data (dict): Nested dictionary {sample_interval: {group_id: group_data}}
        group_number (int): group number
        label (str): label(Such as 'PDB' or 'AF')

    Returns:
        p_values_all (dict): Nested dictionary in the format:
            {
                sample_interval: {
                    'p_levene': np.ndarray(group_number, group_number),
                    'p_ttest': np.ndarray(group_number, group_number)
                },
                ...
            }
    """
    print(f"Performing Levene-tests and T-tests for all intervals in {label}...")

    
    # Define two-dimensional array to save p-values
    p_values_all = {}

    
    # levene-test and t-test for every group and save p-value
    for interval, group_dict in grouped_data.items():
        print(f"\nInterval {interval} ns:")
        p_levenes = np.zeros((group_number, group_number))
        p_ttests = np.zeros((group_number, group_number))
        for i in range(group_number):
            for j in range(group_number):
                data_i = group_dict.get(i, [])
                data_j = group_dict.get(j, [])
                # Ensure there are enough data points
                if len(data_i) == 0 or len(data_j) == 0:
                    print(f"Warning: Empty group at interval {interval}, groups {i}, {j}")
                    p_levenes[i, j] = np.nan
                    p_ttests[i, j] = np.nan
                    continue
                levene_stat, p_levene = stats.levene(data_i, data_j)
                equal_var = True if p_levene > 0.05 else False
                t_stat, p_ttest = stats.ttest_ind(data_i, data_j, equal_var=equal_var)
                
                p_levenes[i, j] = p_levene
                p_ttests[i, j] = p_ttest
                print(f"Group {i+1} vs Group {j+1} — Levene: p={p_levene:.4f}, T-test: p={p_ttest:.4f}")

        # Store the results for that interval
        p_values_all[interval] = {
            'p_levene': p_levenes,
            'p_ttest': p_ttests
        }
    return p_values_all


def plot_internal_p_vs_interval(p_values_all, matrix_key, title, save_path):
    """
    Plot the mean ± standard deviation of the p-value matrix (matrix_key) for group comparisons at different sampling intervals.

    Parameters:
        p_values_all (dict):
            {
                sample_interval: {
                    'p_ttest': np.ndarray(group_number, group_number),
                    'p_levene': np.ndarray(group_number, group_number)
                },
                ...
            }
        matrix_key (str): 'p_ttest' or 'p_levene'
        title (str): Mean ± Std of t-test & levene-test p-values vs Sample Interval
        save_path (str): figure save path
    """
    intervals = sorted(p_values_all.keys())

    mean_vals = []
    std_vals = []

    for interval in intervals:
        matrix = p_values_all[interval][matrix_key]
        # Keep only off-diagonal elements
        upper_tri = matrix[~np.eye(matrix.shape[0], dtype=bool)]
        mean_vals.append(np.mean(upper_tri))
        std_vals.append(np.std(upper_tri))

    # Plot a line graph of the mean ± standard deviation
    plt.figure(figsize=(8, 5))
    plt.errorbar(intervals, mean_vals, yerr=std_vals, fmt='-o', capsize=5, color='blue', ecolor='gray')
    plt.axhline(0.05, color='red', linestyle='--', label='p = 0.05 threshold')
    plt.xlabel('Sample Interval (ns)')
    plt.ylabel(f'Mean ± Std of {matrix_key}')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved lineplot (mean ± std) to {save_path}")


def main():
    args = parse_arguments()

    # Get data
    try:
        print("Reading PDB file...")
        pdb_projection_PC1, pdb_projection_PC2 = read_xvg(args.pdb_file)
        print("Reading AF file...")
        af_projection_PC1, af_projection_PC2 = read_xvg(args.af_file)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)

    # Use sample_data_with_interval to sample the data
    PC1_pdb_af_sampled_data = sample_data_with_interval(
        pdb_projection_PC1, af_projection_PC1, 
        args.min_interval, args.max_interval, args.num_intervals
    )

    # Use group_sampled_data to group the sampled data
    PC1_pdb_group_sampled_data = group_sampled_data(PC1_pdb_af_sampled_data, args.group_number, "PDB")
    PC1_af_group_sampled_data = group_sampled_data(PC1_pdb_af_sampled_data, args.group_number, "AF")

    # Perform tests on the grouped data
    PC1_pdb_internal_p = perform_tests_internal(PC1_pdb_group_sampled_data, args.group_number, "PDB")
    PC1_af_internal_p = perform_tests_internal(PC1_af_group_sampled_data, args.group_number, "AF")

    # Plot figures

    # Plot p-values for PDB internal comparisons
    plot_internal_p_vs_interval(
        PC1_pdb_internal_p,
        'p_ttest',
        f'PDB internal Mean ± Std of t-test p-values vs Sample Interval',
        os.path.join(args.save_dir, f'pdb_internal_t_{args.save_time}_{args.group_number}_PC1.png')
    )
    plot_internal_p_vs_interval(
        PC1_pdb_internal_p,
        'p_levene',
        f'PDB internal Mean ± Std of levene-test p-values vs Sample Interval',
        os.path.join(args.save_dir, f'pdb_internal_levene_{args.save_time}_{args.group_number}_PC1.png')
    )
    # Plot p-values for AF internal comparisons
    plot_internal_p_vs_interval(
        PC1_af_internal_p,
        'p_ttest',
        f'AF Internal Mean ± Std of t-test p-values vs Sample Interval',
        os.path.join(args.save_dir, f'af_internal_t_{args.save_time}_{args.group_number}_PC1.png')
    )
    plot_internal_p_vs_interval(
        PC1_af_internal_p,
        'p_levene',
        f'AF Mean ± Std of levene-test p-values vs Sample Interval',
        os.path.join(args.save_dir, f'af_internal_levene_{args.save_time}_{args.group_number}_PC1.png')
    )
    
        
if __name__ == "__main__":
    main()