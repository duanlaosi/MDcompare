#!/home/duanqi/miniconda3/envs/py36_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
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
                    raise ValueError(f"cannot read {filename}:cannot change '{parts[0]}' or '{parts[1]}' into float \n wrong line:{line.strip()}") from e 
                    


        

    return np.array(x_values), np.array(y_values)

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='PCA Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_file', type=str,
                        help='Path to PDB RMSD of trajectory file')
    
    parser.add_argument('af_file', type=str,
                        help='Path to AF RMSD of trajectory file')
    
    parser.add_argument('plot_RMSD_dir', type=str,
                        help='Directory to save RMSD figure')

    parser.add_argument('save_time', type=str,
                        help='MD length used to name figure')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  Path to PDB RMSD of trajectory file: {args.pdb_file}")
    print(f"  Path to AF RMSD of trajectory file: {args.af_file}")
    print(f"  Directory to save RMSD figure: {args.plot_RMSD_dir}")
    print(f"  MD length used to name figure: {args.save_time}")

    # Get data
    try:
        print("Reading PDB file...")
        pdb_traj_time, pdb_RMSD = read_xvg(args.pdb_file)
        print("Reading AF file...")
        af_traj_time, af_RMSD = read_xvg(args.af_file)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)

    #find max and min in all x(pc1) and y(pc2)
    all_x = np.concatenate([pdb_traj_time, af_traj_time])
    all_y = np.concatenate([pdb_RMSD, af_RMSD])

    #Use downward rounding to get integer minima and upward rounding to get integer maxima
    x_min = int(np.floor(np.min(all_x)))
    x_max = int(np.ceil(np.max(all_x)))
    y_min = 0
    y_max = np.max(all_y)

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(pdb_traj_time, pdb_RMSD, '-', color='grey', linewidth=0.4, label='PDB trajectory RMSD', alpha=0.5)
    plt.plot(af_traj_time, af_RMSD, '-', color='blue', linewidth=0.4, label='AF trajectory RMSD', alpha=0.5)

    #Add title and label
    plt.title('RMSD of PDB trajectory and AF trajectory', fontsize=16)
    plt.xlabel('Time(ns)', fontsize=14)
    plt.ylabel('RMSD(nm)', fontsize=14)

    #plot legend
    plt.legend(fontsize=12)

    #set x axis ticks
    num_ticks = 11
    step = step=max(1, int(np.round((x_max - x_min) / (num_ticks - 1))))
    xticks = np.arange(x_min, x_max, step)
    plt.xticks(xticks)

    #Setting x axis range
    plt.xlim(x_min - 1, x_max + 1)

    #set y axis ticks
    step=0.2
    # num_steps = int(np.ceil((y_max - y_min) / step))
    # y_max_adjusted = y_min + num_steps * step

    yticks = np.arange(y_min, y_max + step, step)
    plt.yticks(yticks)

    #Setting y axis range
    plt.ylim(y_min, y_max + step)

    #Setting plot path
    plot_path = os.path.join(args.plot_RMSD_dir, f'RMSD_{args.save_time}.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Saved RMSD plot to {plot_path}")

if __name__ == "__main__":
    main()