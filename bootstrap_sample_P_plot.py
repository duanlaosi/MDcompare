#!/home/duanqi/miniconda3/envs/py36_env/bin/python
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

def bootstrap_ttest_pval(data1, data2, n_boot=1000):
    """
    Perform bootstrap resampling to get t test p-value distribution
    """
    pvals = []
    for _ in range(n_boot):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        _, p = stats.ttest_ind(sample1, sample2, equal_var=False)
        pvals.append(p)
    
    pvals = np.array(pvals)
    return np.median(pvals), np.percentile(pvals, [2.5, 97.5])

def bootstrap_levene_pval(data1, data2, n_boot=1000):
    """
    Perform bootstrap resampling to get Levene test p-value distribution
    """
    pvals = []
    for _ in range(n_boot):
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        _, p = stats.levene(sample1, sample2)
        pvals.append(p)
    
    pvals = np.array(pvals)

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
            'ci_upper': ttest_ci[1],
            'pvals': ttest_pvals
        },
        'levene': {
            'median': levene_median,
            'mean': levene_mean,
            'ci_lower': levene_ci[0],
            'ci_upper': levene_ci[1],
            'pvals': levene_pvals
        }
    }

def plot_p_distribution(pvals, test_name, interval, save_dir, save_time):
    """
    Plot distribution of bootstrap p-values
    """
    plt.figure(figsize=(8, 5))
    plt.hist(pvals, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
    plt.axvline(0.05, color='red', linestyle='--', label='p = 0.05')
    plt.xlabel('p-value')
    plt.ylabel('Frequency')
    plt.title(f'{test_name} p-value Distribution (Interval = {interval:.1f} ns)')
    plt.legend()
    plt.tight_layout()
    fname = os.path.join(save_dir, f'pdist_{test_name.lower()}_{interval:.1f}ns_{save_time}.png')
    plt.savefig(fname, dpi=300)
    plt.close()


def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  PDB file: {args.pdb_file}")
    print(f"  AF file: {args.af_file}")
    print(f"  Min interval: {args.min_interval}")
    print(f"  Max interval: {args.max_interval}")
    print(f"  Number of intervals: {args.num_intervals}")
    print(f"  Save directory: {args.save_dir}")
    print(f"  Bootstrap iterations: {args.n_boot}")
    print(f"  Confidence level: {args.confidence_level}%")
    
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





    #generate a A list of equally spaced sampling intervals(useing num_intervals)
    # Generate sampling intervals
    min_int = int(args.min_interval * 10)
    max_int = int(args.max_interval * 10)
    step = (max_int - min_int) // (args.num_intervals - 1)
    int_intervals = list(range(min_int, max_int + 1, step))
    sample_intervals = [val / 10 for val in int_intervals]
    step_sizes = [max(1, int_val) for int_val in int_intervals]

    # Initialize storage for results
    original_results = {'levene': [], 'ttest': []}
    bootstrap_results = {'levene': [], 'ttest': []}

    print("\nCalculating statistics with bootstrap...")
   
   # Loop through each sampling interval
    for i, sample_interval in enumerate(sample_intervals):
        step = step_sizes[i]
        print(f"Processing interval {sample_interval:.1f}ns (step {step})...")
        
        # Sample the data
        pdb_sample_PC1 = pdb_projection_PC1[::step]
        af_sample_PC1 = af_projection_PC1[::step]
        
        # Original (non-bootstrap) analysis
        stat_levene, p_levene = stats.levene(af_sample_PC1, pdb_sample_PC1)
        equal_var = True if p_levene > 0.05 else False 
        t_stat, p_ttest = stats.ttest_ind(af_sample_PC1, pdb_sample_PC1, equal_var=equal_var)
        
        original_results['levene'].append(p_levene)
        original_results['ttest'].append(p_ttest)
        
        # Bootstrap analysis
        boot_stats = bootstrap_statistics(af_sample_PC1, pdb_sample_PC1, 
                                        args.n_boot, args.confidence_level)
        
        # Plot p-value distributions for this interval
        plot_p_distribution(boot_stats['ttest']['pvals'], 't-test', sample_interval, args.save_dir, args.save_time)
        plot_p_distribution(boot_stats['levene']['pvals'], 'Levene', sample_interval, args.save_dir, args.save_time)



        bootstrap_results['levene'].append(boot_stats['levene'])
        bootstrap_results['ttest'].append(boot_stats['ttest'])
        
        print(f"  Original: Levene p={p_levene:.4f}, t-test p={p_ttest:.4f}")
        print(f"  Bootstrap: Levene median={boot_stats['levene']['median']:.4f}, t-test median={boot_stats['ttest']['median']:.4f}")

    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    print("\nGenerating plots...")
    
    # Extract bootstrap results for plotting
    levene_medians = [result['median'] for result in bootstrap_results['levene']]
    levene_means = [result['mean'] for result in bootstrap_results['levene']]
    levene_ci_lower = [result['ci_lower'] for result in bootstrap_results['levene']]
    levene_ci_upper = [result['ci_upper'] for result in bootstrap_results['levene']]
    
    ttest_medians = [result['median'] for result in bootstrap_results['ttest']]
    ttest_means = [result['mean'] for result in bootstrap_results['ttest']]
    ttest_ci_lower = [result['ci_lower'] for result in bootstrap_results['ttest']]
    ttest_ci_upper = [result['ci_upper'] for result in bootstrap_results['ttest']]

    # Plot Levene test results
    plt.figure(figsize=(12, 8))

    # Plot confidence interval band
    plt.fill_between(sample_intervals, levene_ci_lower, levene_ci_upper, 
                    alpha=0.3, color='lightblue', label=f'{args.confidence_level}% CI')
    
    # Plot bootstrap median/mean
    plt.plot(sample_intervals, levene_medians, 'b-', linewidth=2, 
            label='Bootstrap Median', marker='s', markersize=6)
    plt.plot(sample_intervals, levene_means, 'b--', linewidth=2, 
            label='Bootstrap Mean', marker='^', markersize=6)
    
    # Plot original data points
    plt.plot(sample_intervals, original_results['levene'], 'ro', 
            markersize=8, label='Original Data', markerfacecolor='red', 
            markeredgecolor='darkred', markeredgewidth=2)

    # Add significance threshold
    plt.axhline(0.05, color='red', linestyle='--', alpha=0.7, label='p = 0.05 threshold')


    plt.xlabel('Sample Interval (ns)', fontsize=12)
    plt.ylabel('p-value (Levene Test)', fontsize=12)
    plt.title(f'Levene Test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)', fontsize=14)
    plt.xticks(sample_intervals)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    levene_plot_path = os.path.join(args.save_dir, f'bootstrap_levene_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(levene_plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved bootstrap Levene test plot to {levene_plot_path}")
    
    # Plot t-test results
    plt.figure(figsize=(12, 8))
    
    # Plot confidence interval band
    plt.fill_between(sample_intervals, ttest_ci_lower, ttest_ci_upper, 
                    alpha=0.3, color='lightgreen', label=f'{args.confidence_level}% CI')
    
    # Plot bootstrap median/mean
    plt.plot(sample_intervals, ttest_medians, 'g-', linewidth=2, 
            label='Bootstrap Median', marker='s', markersize=6)
    plt.plot(sample_intervals, ttest_means, 'g--', linewidth=2, 
            label='Bootstrap Mean', marker='^', markersize=6)
    
    # Plot original data points
    plt.plot(sample_intervals, original_results['ttest'], 'ro', 
            markersize=8, label='Original Data', markerfacecolor='red', 
            markeredgecolor='darkred', markeredgewidth=2)

    # Add significance threshold
    plt.axhline(0.05, color='red', linestyle='--', alpha=0.7, label='p = 0.05 threshold')

    plt.xlabel('Sample Interval (ns)', fontsize=12)
    plt.ylabel('p-value (t-test)', fontsize=12)
    plt.title(f't-test p-values vs Sample Interval\n(Bootstrap: {args.n_boot} iterations)', fontsize=14)
    plt.xticks(sample_intervals)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    ttest_plot_path = os.path.join(args.save_dir, f'bootstrap_ttest_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(ttest_plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved bootstrap t-test plot to {ttest_plot_path}")
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()


