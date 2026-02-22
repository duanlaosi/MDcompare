#!/home/duanqi/miniconda3/envs/py38_env/bin/python
# Two PC ANOVA (Hotelling's T-squared test) plot
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import os
import sys
import argparse
import pandas as pd
import pingouin as pg
from pingouin import multivariate_ttest
# import inspect
# print(inspect.signature(pg.box_m))
# print(pg.__version__)
# print(pg.box_m.__module__)
# print(pg.box_m.__code__.co_filename)


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
    
    return parser.parse_args()

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




    #plot test(variance test)'s p-value of sample interval
    p_boxm = []
    p_manova = []
    #generate a A list of equally spaced sampling intervals(useing num_intervals)
    min_int = int(args.min_interval * 10)
    max_int = int(args.max_interval * 10)
    step = (max_int - min_int) // (args.num_intervals - 1)
    int_intervals = list(range(min_int, max_int + 1, step))
    sample_intervals = [val / 10 for val in int_intervals]
    step_sizes = [max(1, int_val) for int_val in int_intervals]

    print("\nCalculating statistics...")
    # Loop through stepsize from 1 to 10(0.1ns sample interval to 1ns sample interval correspond to trajectory frame 1 to 10)
    for i, sample_interval in enumerate(sample_intervals):
        step = step_sizes[i]
        # Calculate the sampling for PDB data
        pdb_sample_PC1 = pdb_projection_PC1[::step]
        pdb_sample_PC2 = pdb_projection_PC2[::step]
        
        # Calculate the sampling for AF data
        af_sample_PC1 = af_projection_PC1[::step]
        af_sample_PC2 = af_projection_PC2[::step]

        PC1 = np.concatenate([pdb_sample_PC1, af_sample_PC1])
        PC2 = np.concatenate([pdb_sample_PC2, af_sample_PC2])

        groups = ['PDB'] * len(pdb_sample_PC1) + ['AF'] * len(af_sample_PC1)

        df = pd.DataFrame({
            'PC1': PC1,
            'PC2': PC2,
            'group': groups
        })

        box_result = pg.box_m(data=df, dvs=['PC1', 'PC2'], group='group')
        p_box = box_result['pval'].iloc[0]
        p_boxm.append(p_box)

        if p_box > 0.05:
            method = 'Wilks'
        else:
            method = 'Pillai'

        # Perform Hotelling's T-squared test (MANOVA)
        X = np.column_stack([af_sample_PC1, af_sample_PC2])
        Y = np.column_stack([pdb_sample_PC1, pdb_sample_PC2])
        test_result = multivariate_ttest(X=X, Y=Y, paired=False)
        p_val = test_result['pval'].iloc[0]
        p_manova.append(p_val)

        print(f"Sample interval: {sample_interval:.1f}, Step size: {step}, Box's M p-value: {p_box:.4f}, MANOVA ({method}) p-value: {p_val:.4f}")

    os.makedirs(args.save_dir, exist_ok=True)
    
    print("\nGenerating plots...")
    # Create a plot for Box's M test p-values 
    plt.figure(figsize=(10, 6))
    plt.plot(sample_intervals, p_boxm, 'o-', markeredgecolor='blue', markerfacecolor='none', markersize=8, clip_on=False)
    plt.axhline(0.05, color='red', linestyle='--', label='Threshold 0.05')
    plt.xlabel('Sample Interval', fontsize=12)
    plt.ylabel("p-value (Box's M Test)", fontsize=12)
    plt.title("Box's M p-values vs Sample Interval", fontsize=14)
    plt.xticks(sample_intervals)
    for i, txt in enumerate(p_boxm):
        plt.annotate(f"{txt:.4f}", (sample_intervals[i], p_boxm[i]), textcoords="offset points", xytext=(0,10), ha='center')
    plt.legend()
    plt.tight_layout()
    box_plot_path = os.path.join(args.save_dir, f'p_boxm_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(box_plot_path)
    print(f"Saved Box's M test plot to {box_plot_path}")

    # Create a plot for MANOVA p-values
    plt.figure(figsize=(10, 6))
    plt.plot(sample_intervals, p_manova, 'o-', markeredgecolor='green', markerfacecolor='none', markersize=8, clip_on=False)
    plt.axhline(0.05, color='red', linestyle='--', label='Threshold 0.05')
    plt.xlabel('Sample Interval', fontsize=12)
    plt.ylabel('p-value (MANOVA)', fontsize=12)
    plt.title('MANOVA p-values vs Sample Interval', fontsize=14)
    plt.xticks(sample_intervals)
    for i, txt in enumerate(p_manova):
        plt.annotate(f"{txt:.4f}", (sample_intervals[i], p_manova[i]), textcoords="offset points", xytext=(0,10), ha='center')
    plt.legend()
    plt.tight_layout()
    manova_plot_path = os.path.join(args.save_dir, f'p_manova_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(manova_plot_path)
    print(f"Saved MANOVA test plot to {manova_plot_path}")
if __name__ == "__main__":
    main()


