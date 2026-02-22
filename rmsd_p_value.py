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
    print(f"[DEBUG] Entering read_xvg with file: {filename}")
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
    parser = argparse.ArgumentParser(description='Analysis Script')
    
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
        sys.stdout.flush()
        pdb_time, pdb_RMSD = read_xvg(args.pdb_file)
        print(f"Reading PDB file...{args.pdb_file}")
        print(f"PDB RMSD values (first 5): {pdb_RMSD[:5]}")
        print(f"Total PDB data points: {len(pdb_RMSD)}")
        print(f"[DEBUG] args.pdb_file = {args.pdb_file}")
        sys.stdout.flush()

        print("Reading AF file...")
        sys.stdout.flush()
        af_time, af_RMSD = read_xvg(args.af_file)
        print(f"Reading AF file...{args.af_file}")
        print(f"AF RMSD values (first 5): {af_RMSD[:5]}")
        print(f"Total AF data points: {len(af_RMSD)}")
        print(f"[DEBUG] args.af_file = {args.af_file}")
    except Exception as e:
        import traceback
        print(f"[ERROR] Exception occurred while reading XVG files:\n{e}")
        traceback.print_exc()
        sys.exit(1)




    #plot t-test and levene test(variance test)'s p-value from 0.1ns sample interval to 10ns sample interval
    p_levenes = []
    p_ttests = []
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
        pdb_sample_time = pdb_time[::step]
        pdb_sample_RMSD = pdb_RMSD[::step]
        
        # Calculate the sampling for AF data
        af_sample_time = af_time[::step]
        af_sample_RMSD = af_RMSD[::step]

        stat_levene, p_levene = stats.levene(af_sample_RMSD, pdb_sample_RMSD)
        equal_var = True if p_levene > 0.05 else False 
        t_stat, p_ttest = stats.ttest_ind(af_sample_RMSD, pdb_sample_RMSD, equal_var=equal_var)
        
        # Store the p-values
        p_levenes.append(p_levene)
        p_ttests.append(p_ttest)

        print(f"Sample interval: {sample_interval:.1f}, Step size: {step}")

    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    print("\nGenerating plots...")
    # Create a plot for p_levenes
    plt.figure(figsize=(10, 6))
    plt.plot(sample_intervals, p_levenes, 'o-', markeredgecolor='blue', markerfacecolor='none', markersize=8, clip_on=False)
    plt.xlabel('Sample Interval', fontsize=12)
    plt.ylabel('p-value (Levene Test)', fontsize=12)
    plt.title('Levene Test p-values vs Sample Interval', fontsize=14)
    plt.xticks(sample_intervals)
    for i, txt in enumerate(p_levenes):
        plt.annotate(f"{txt:.4f}", (sample_intervals[i], p_levenes[i]), 
                    textcoords="offset points", xytext=(0,10), ha='center')
    plt.tight_layout()
    levene_plot_path = os.path.join(args.save_dir, f'p_levenes_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(levene_plot_path)
    print(f"Saved Levene test plot to {levene_plot_path}")

    # Create a plot for p_ttests
    plt.figure(figsize=(10, 6))
    plt.plot(sample_intervals, p_ttests, 'o-', markeredgecolor='blue', markerfacecolor='none', markersize=8, clip_on=False)
    plt.xlabel('Sample Interval', fontsize=12)
    plt.ylabel('p-value (t-test)', fontsize=12)
    plt.title('t-test p-values vs Sample Interval', fontsize=14)
    plt.xticks(sample_intervals)
    for i, txt in enumerate(p_ttests):
        plt.annotate(f"{txt:.4f}", (sample_intervals[i], p_ttests[i]), 
                    textcoords="offset points", xytext=(0,10), ha='center')
    plt.tight_layout()
    ttest_plot_path = os.path.join(args.save_dir, f'p_ttests_{args.min_interval}_{args.max_interval}_{args.save_time}.png')
    plt.savefig(ttest_plot_path)
    print(f"Saved t-test plot to {ttest_plot_path}")
if __name__ == "__main__":
    main()


