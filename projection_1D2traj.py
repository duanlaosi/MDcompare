#!/home/duanqi/miniconda3/envs/py36_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

def read_xvg_gX(filename):
    """
    Read .xvg files @ with gX groups
    Returns each set of data as a list, where each element is a tuple (x_values, y_values)
    """
    data_groups = []
    x_values = []
    y_values = []
    in_group = False
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            # skip title and annotations and empty line
            if line.startswith('&'):
                continue
            if not line:
                continue
            
            if line.startswith('@ with g'):
                # If there is already a group being read, save it
                if in_group and x_values:
                    data_groups.append((np.array(x_values), np.array(y_values)))
                    x_values = []
                    y_values = []
                in_group = True
                continue
            
            if line.startswith('@'):
                continue

            # split data lines
            parts = line.split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    x_values.append(x)
                    y_values.append(y)
                except ValueError as e:
                    raise ValueError(f"cannot read {filename}:cannot change '{parts[0]}' or '{parts[1]}' into float \n wrong line:{line}") from e 
                    
    if in_group and x_values:
        data_groups.append((np.array(x_values), np.array(y_values)))
    
    return data_groups

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='PCA Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_file', type=str,
                        help='Path to PDB projection of trajectory file')
    
    parser.add_argument('af_file', type=str,
                        help='Path to AF projection of trajectory file')
    
    parser.add_argument('plot_PCA_dir', type=str,
                        help='Directory to save projection figure')

    parser.add_argument('save_time', type=str,
                        help='MD length used to name figure')
    
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  Path to PDB projection of trajectory file: {args.pdb_file}")
    print(f"  Path to AF projection of trajectory file: {args.af_file}")
    print(f"  Directory to save projection figure: {args.plot_PCA_dir}")
    print(f"  MD length used to name figure: {args.save_time}")

    # Get data
    try:
        print("Reading PDB file...")
        groups_pdb = read_xvg_gX(args.pdb_file)
        pdb_traj_time_PC1, pdb_proj_PC1 = groups_pdb[0]
        print("Reading AF file...")
        groups_af = read_xvg_gX(args.af_file)
        af_traj_time_PC1, af_proj_PC1 = groups_af[0]
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)

    # find max and min in all x(pc1) and y(pc2)
    all_x = np.concatenate([pdb_traj_time_PC1, af_traj_time_PC1])
    all_y = np.concatenate([pdb_proj_PC1, af_proj_PC1])

    # Use downward rounding to get integer minima and upward rounding to get integer maxima
    x_min = int(np.floor(np.min(all_x)))
    x_max = int(np.ceil(np.max(all_x)))
    y_min = np.floor(np.min(all_y))
    y_max = np.ceil(np.max(all_y))

    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(pdb_traj_time_PC1, pdb_proj_PC1, '-', color='black', linewidth=0.4, label='PDB trajectory Projection', alpha=0.5)
    plt.plot(af_traj_time_PC1, af_proj_PC1, '-', color='blue', linewidth=0.4, label='AF trajectory Projection', alpha=0.5)

    # Add title and label
    plt.title('PC1 projection of PDB trajectory and AF trajectory', fontsize=16)
    plt.xlabel('Time(ns)', fontsize=14)
    plt.ylabel('Projection on PC1(nm)', fontsize=14)

    # plot legend
    plt.legend(fontsize=12)

    # set x axis ticks
    num_ticks = 11
    step = step=max(1, int(np.round((x_max - x_min) / (num_ticks - 1))))
    xticks = np.arange(x_min, x_max + 1, step)
    plt.xticks(xticks)

    # Setting x axis range
    plt.xlim(x_min - 1, x_max + 1)

    # set y axis ticks
    num_ticks = 10
    step = step=max(1, int(np.round((y_max - y_min) / (num_ticks - 1))))
    num_steps = (y_max - y_min) // step
    y_max_adjusted = y_min + num_steps * step
    if y_max_adjusted < y_max:
        y_max_adjusted += step
    yticks = np.arange(y_min, y_max_adjusted + step, step)
    plt.yticks(yticks)

    # Setting y axis range
    plt.ylim(y_min, y_max_adjusted)

    plt.tight_layout()

    # Setting plot path
    plot_path = os.path.join(args.plot_PCA_dir, f'PC1_proj_{args.save_time}.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Saved projection plot to {plot_path}")

if __name__ == "__main__":
    main()