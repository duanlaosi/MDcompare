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
    
    # Add optional bootstrap parameters
    parser.add_argument('--n_boot', type=int, default=1000,
                        help='Number of bootstrap iterations (default: 1000)')
    
    parser.add_argument('--confidence_level', type=float, default=95,
                        help='Confidence level for intervals (default: 95)')
    
    return parser.parse_args()

def bootstrap_statistics(data1, data2, n_boot=1000, confidence_level=95):
    """
    Calculate bootstrap statistics for both t-test and Levene test
    """
    alpha = (100 - confidence_level) / 2
    percentiles = [alpha, 100 - alpha]
    
    # Bootstrap t-test
    ttest_pvals = []
    levene_pvals = []
    
    for _ in range(n_boot):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        
        # t-test
        _, p_ttest = stats.ttest_ind(sample1, sample2, equal_var=False)
        ttest_pvals.append(p_ttest)
        
        # Levene test
        _, p_levene = stats.levene(sample1, sample2)
        levene_pvals.append(p_levene)
    
    ttest_pvals = np.array(ttest_pvals)
    levene_pvals = np.array(levene_pvals)
    
    # Calculate statistics
    ttest_median = np.median(ttest_pvals)
    ttest_mean = np.mean(ttest_pvals)
    ttest_ci = np.percentile(ttest_pvals, percentiles)
    
    levene_median = np.median(levene_pvals)
    levene_mean = np.mean(levene_pvals)
    levene_ci = np.percentile(levene_pvals, percentiles)
    
    return {
        'ttest': {
            'median': ttest_median,
            'mean': ttest_mean,
            'ci_lower': ttest_ci[0],
            'ci_upper': ttest_ci[1]
        },
        'levene': {
            'median': levene_median,
            'mean': levene_mean,
            'ci_lower': levene_ci[0],
            'ci_upper': levene_ci[1]
        }
    }


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

def perform_bootstrap_tests_internal(grouped_data, group_number, label, n_boot=1000, confidence_level=95):
    """
    Perform bootstrap resampling for Levene's test and t-test on grouped data at multiple sampling intervals
    
    Parameters:
        grouped_data (dict): Nested dictionary {sample_interval: {group_id: group_data}}
        group_number (int): group number
        label (str): label (Such as 'PDB' or 'AF')
        n_boot (int): number of bootstrap iterations
        confidence_level (float): confidence level for intervals
    
    Returns:
        bootstrap_results (dict): Nested dictionary with bootstrap statistics
    """
    print(f"Performing bootstrap analysis for {label} with {n_boot} iterations...")
    
    bootstrap_results = {}
    
    for interval, group_dict in grouped_data.items():
        print(f"\nBootstrap analysis for interval {interval} ns:")
        
        # Collect all pairwise comparisons for this interval
        interval_ttest_stats = []
        interval_levene_stats = []
        
        for i in range(group_number):
            for j in range(i+1, group_number):  # Only upper triangle to avoid duplicates
                data_i = group_dict.get(i, [])
                data_j = group_dict.get(j, [])
                
                if len(data_i) == 0 or len(data_j) == 0:
                    print(f"Warning: Empty group at interval {interval}, groups {i}, {j}")
                    continue
                
                # Perform bootstrap analysis for this pair
                bootstrap_stats = bootstrap_statistics(data_i, data_j, n_boot, confidence_level)
                interval_ttest_stats.append(bootstrap_stats['ttest'])
                interval_levene_stats.append(bootstrap_stats['levene'])
        
        # Calculate overall statistics for this interval
        if interval_ttest_stats:
            # Extract medians and means from all pairwise comparisons
            ttest_medians = [stat['median'] for stat in interval_ttest_stats]
            ttest_means = [stat['mean'] for stat in interval_ttest_stats]
            ttest_ci_lowers = [stat['ci_lower'] for stat in interval_ttest_stats]
            ttest_ci_uppers = [stat['ci_upper'] for stat in interval_ttest_stats]
            
            levene_medians = [stat['median'] for stat in interval_levene_stats]
            levene_means = [stat['mean'] for stat in interval_levene_stats]
            levene_ci_lowers = [stat['ci_lower'] for stat in interval_levene_stats]
            levene_ci_uppers = [stat['ci_upper'] for stat in interval_levene_stats]
            
            # Store results for this interval
            bootstrap_results[interval] = {
                'ttest': {
                    'overall_median': np.median(ttest_medians),
                    'overall_mean': np.mean(ttest_means),
                    'ci_lower': np.mean(ttest_ci_lowers),
                    'ci_upper': np.mean(ttest_ci_uppers),
                    'original_data': np.mean(ttest_medians)  # Use median of bootstrap medians as original
                },
                'levene': {
                    'overall_median': np.median(levene_medians),
                    'overall_mean': np.mean(levene_means),
                    'ci_lower': np.mean(levene_ci_lowers),
                    'ci_upper': np.mean(levene_ci_uppers),
                    'original_data': np.mean(levene_medians)  # Use median of bootstrap medians as original
                }
            }
            
            print(f"  T-test - Bootstrap median: {bootstrap_results[interval]['ttest']['overall_median']:.4f}, "
                  f"95% CI: [{bootstrap_results[interval]['ttest']['ci_lower']:.4f}, "
                  f"{bootstrap_results[interval]['ttest']['ci_upper']:.4f}]")
            print(f"  Levene - Bootstrap median: {bootstrap_results[interval]['levene']['overall_median']:.4f}, "
                  f"95% CI: [{bootstrap_results[interval]['levene']['ci_lower']:.4f}, "
                  f"{bootstrap_results[interval]['levene']['ci_upper']:.4f}]")
    
    return bootstrap_results

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

    
    # Define dictionary to save p-values
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


def plot_bootstrap_results(bootstrap_results, test_type, title, save_path):
    """
    Plot bootstrap results with confidence intervals
    
    Parameters:
        bootstrap_results (dict): Bootstrap analysis results
        test_type (str): 'ttest' or 'levene'
        title (str): Plot title
        save_path (str): Save path for the plot
    """
    intervals = sorted(bootstrap_results.keys())
    
    medians = []
    means = []
    ci_lowers = []
    ci_uppers = []
    original_data = []
    
    for interval in intervals:
        stats = bootstrap_results[interval][test_type]
        medians.append(stats['overall_median'])
        means.append(stats['overall_mean'])
        ci_lowers.append(stats['ci_lower'])
        ci_uppers.append(stats['ci_upper'])
        original_data.append(stats['original_data'])
    
    medians = np.array(medians)
    means = np.array(means)
    ci_lowers = np.array(ci_lowers)
    ci_uppers = np.array(ci_uppers)
    original_data = np.array(original_data)
    
    plt.figure(figsize=(10, 6))
    
    # Plot confidence interval as filled area
    plt.fill_between(intervals, ci_lowers, ci_uppers, alpha=0.3, color='lightblue', label='95% CI')
    
    # Plot bootstrap median and mean
    plt.plot(intervals, medians, 'b-', marker='s', linewidth=2, label='Bootstrap Median')
    plt.plot(intervals, means, 'b--', marker='^', linewidth=2, label='Bootstrap Mean')
    
    # Plot original data points
    plt.plot(intervals, original_data, 'ro', markersize=8, label='Original Data')
    
    # Add significance threshold
    plt.axhline(0.05, color='red', linestyle='--', alpha=0.7, label='p = 0.05 threshold')
    
    plt.xlabel('Sample Interval (ns)')
    plt.ylabel(f'p-value ({test_type.replace("_", " ").title()})')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved bootstrap plot to {save_path}")



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

    # Perform original tests on the grouped data
    PC1_pdb_internal_p = perform_tests_internal(PC1_pdb_group_sampled_data, args.group_number, "PDB")
    PC1_af_internal_p = perform_tests_internal(PC1_af_group_sampled_data, args.group_number, "AF")

    # Perform bootstrap analysis
    print(f"\nStarting bootstrap analysis with {args.n_boot} iterations...")
    PC1_pdb_bootstrap = perform_bootstrap_tests_internal(
        PC1_pdb_group_sampled_data, args.group_number, "PDB", 
        args.n_boot, args.confidence_level
    )
    PC1_af_bootstrap = perform_bootstrap_tests_internal(
        PC1_af_group_sampled_data, args.group_number, "AF", 
        args.n_boot, args.confidence_level
    )

    # Plot bootstrap results
    plot_bootstrap_results(
        PC1_pdb_bootstrap,
        'ttest',
        f'PDB Internal T-test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)',
        os.path.join(args.save_dir, f'pdb_internal_t_bootstrap_{args.save_time}_{args.group_number}_PC1.png')
    )
    plot_bootstrap_results(
        PC1_pdb_bootstrap,
        'levene',
        f'PDB Internal Levene Test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)',
        os.path.join(args.save_dir, f'pdb_internal_levene_bootstrap_{args.save_time}_{args.group_number}_PC1.png')
    )
    plot_bootstrap_results(
        PC1_af_bootstrap,
        'ttest',
        f'AF Internal T-test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)',
        os.path.join(args.save_dir, f'af_internal_t_bootstrap_{args.save_time}_{args.group_number}_PC1.png')
    )
    plot_bootstrap_results(
        PC1_af_bootstrap,
        'levene',
        f'AF Internal Levene Test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)',
        os.path.join(args.save_dir, f'af_internal_levene_bootstrap_{args.save_time}_{args.group_number}_PC1.png')
    )
    
    print(f"\nBootstrap analysis completed. All plots saved to {args.save_dir}")
        
if __name__ == "__main__":
    main()