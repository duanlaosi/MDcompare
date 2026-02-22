#!/home/duanqi/miniconda3/envs/py38_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
import sys
import argparse
from scipy.ndimage import gaussian_filter

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
    parser = argparse.ArgumentParser(description='Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_temperature_traj', type=str,
                        help='Path to PDB temperature of trajectory file')

    parser.add_argument('pdb_proj_traj', type=str,
                        help='Path to PDB projection of trajectory file')
    
    parser.add_argument('pdb_proj_nvt', type=str,
                    help='Path to PDB projection of nvt file')
    
    parser.add_argument('pdb_proj_em', type=str,
                help='Path to PDB projection of em file')
    
    parser.add_argument('af_temperature_traj', type=str,
                    help='Path to AF temperature of trajectory file')
    
    parser.add_argument('af_proj_traj', type=str,
                        help='Path to AF projection of trajectory file')

    parser.add_argument('af_proj_nvt', type=str,
                        help='Path to AF projection of nvt file')

    parser.add_argument('af_proj_em', type=str,
                help='Path to AF projection of em file')   
    
    parser.add_argument('save_dir', type=str,
                        help='Directory to save the plots')
    
    parser.add_argument('save_time', type=str,
                        help='Time to save the plots')
    
    return parser.parse_args()

# Function to calculate and plot free energy surface
def plot_free_energy_surface(pc1_data, pc2_data, T, k_B, title, save_dir, pdb_nvt_projection_PC1, pdb_nvt_projection_PC2, pdb_em_projection_PC1, pdb_em_projection_PC2, af_nvt_projection_PC1, af_nvt_projection_PC2, af_em_projection_PC1, af_em_projection_PC2, save_time, add_structures=True):
    """
    Calculate and plot free energy surface
    """
    #Use downward rounding to get integer minima and upward rounding to get integer maxima
    x_min = np.floor(np.min(pc1_data))
    x_max = np.ceil(np.max(pc1_data))
    y_min = np.floor(np.min(pc2_data))
    y_max = np.ceil(np.max(pc2_data))


    # 计算2D直方图（概率密度）
    nbins = 100
    H, xedges, yedges = np.histogram2d(pc1_data, pc2_data, bins=nbins, 
                                      range=[[x_min, x_max], [y_min, y_max]], 
                                      density=True)
    
    # 避免 log(0)
    H[H == 0] = 1e-12
    H = gaussian_filter(H, sigma=1)

    # 计算ΔG (自由能)
    F = -k_B * T * np.log(H)
    F = F - np.min(F)  # 自由能归一化
    
    # 绘图
    X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
    plt.figure(figsize=(10, 8))
    cmap = plt.get_cmap("viridis")
    contour = plt.contourf(X, Y, F.T, levels=50, cmap=cmap)
    cbar = plt.colorbar(contour)
    cbar.set_label('Free Energy (kJ/mol)')

    # Add structure markers if requested
    if add_structures:
        # Plot PDB structures (black markers)
        plt.plot(pdb_nvt_projection_PC1, pdb_nvt_projection_PC2, 'o', 
                color='black', markersize=8, alpha=0.8, markeredgecolor='white', 
                markeredgewidth=1, label='PDB NVT')
        plt.plot(pdb_em_projection_PC1, pdb_em_projection_PC2, 's', 
                color='black', markersize=8, alpha=0.8, markeredgecolor='white',
                markeredgewidth=1, label='PDB EM')
        
        # Plot AF structures (blue markers)
        plt.plot(af_nvt_projection_PC1, af_nvt_projection_PC2, 'o', 
                color='blue', markersize=8, alpha=0.8, markeredgecolor='white',
                markeredgewidth=1, label='AF NVT')
        plt.plot(af_em_projection_PC1, af_em_projection_PC2, 's', 
                color='blue', markersize=8, alpha=0.8, markeredgecolor='white',
                markeredgewidth=1, label='AF EM')
        
        # Smart positioning to avoid overlap
        offset_x = 0.1 * (x_max - x_min)
        offset_y = 0.1 * (y_max - y_min)
        # Add annotations with better positioning
        annotations = [
            (pdb_nvt_projection_PC1[0], pdb_nvt_projection_PC2[0], 'PDB after equilibration', 'black', offset_x, offset_y),
            (pdb_em_projection_PC1[0], pdb_em_projection_PC2[0], 'PDB after energy minimization', 'black', -offset_x, offset_y),
            (af_nvt_projection_PC1[0], af_nvt_projection_PC2[0], 'AF after equilibration', 'blue', offset_x, -offset_y),
            (af_em_projection_PC1[0], af_em_projection_PC2[0], 'AF after energy minimization', 'blue', -offset_x, -offset_y)
        ]
        
        for x, y, text, color, dx, dy in annotations:
            
            plt.annotate(text,
                        xy=(x, y),               
                        xytext=(x + dx, y + dy),   
                        arrowprops=dict(arrowstyle="->", color=color, alpha=0.7),
                        fontsize=6, ha='left',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
        
        # Add legend
        plt.legend(loc='upper right', framealpha=0.9)

        
        
    plt.xlabel('PC1', fontsize=14)
    plt.ylabel('PC2', fontsize=14)
    plt.title(f'Free Energy Surface on {title}', fontsize=16)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.tight_layout()
    
    #set x axis ticks
    num_ticks = 10
    step=max(1, int(np.round((x_max - x_min) / (num_ticks - 1))))
    num_steps = (x_max - x_min) // step
    x_max_adjusted = x_min + num_steps * step
    if x_max_adjusted < x_max:
        x_max_adjusted += step

    xticks = np.arange(x_min, x_max_adjusted + step, step)
    plt.xticks(xticks)


    #set y axis ticks
    num_ticks = 10
    step=max(1, int(np.round((y_max - y_min) / (num_ticks - 1))))
    num_steps = (y_max - y_min) // step
    y_max_adjusted = y_min + num_steps * step
    if y_max_adjusted < y_max:
        y_max_adjusted += step

    yticks = np.arange(y_min, y_max_adjusted + step, step)
    plt.yticks(yticks)

    # Save plot
    save_path = os.path.join(save_dir, f"2traj_add_str{add_structures}_FEP_{title}_{save_time}.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.savefig(save_path.replace(".png", ".pdf"), dpi=300, bbox_inches='tight')

    
    print(f"Free energy surface saved: {save_path}")
    print(f"Max free energy: {np.max(F):.2f} kJ/mol")
    
    return F, X, Y

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  PDB temperature of trajectory file: {args.pdb_temperature_traj}")
    print(f"  PDB projection of trajectory file: {args.pdb_proj_traj}")
    print(f"  PDB projection of nvt file: {args.pdb_proj_nvt}")
    print(f"  PDB projection of em file: {args.pdb_proj_em}")
    print(f"  AF temperature of trajectory file: {args.af_temperature_traj}")
    print(f"  AF projection of trajectory file: {args.af_proj_traj}")
    print(f"  AF projection of nvt file: {args.af_proj_nvt}")
    print(f"  AF projection of em file: {args.af_proj_em}")
    print(f"  Save directory: {args.save_dir}")
    
    # Check if the files exist

    if not os.path.exists(args.pdb_temperature_traj):
        print(f"Error: PDB temperature of trajectory file '{args.pdb_temperature_traj}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_traj):
        print(f"Error: PDB projection of trajectory file '{args.pdb_proj_traj}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_nvt):
        print(f"Error: PDB projection of nvt file '{args.pdb_proj_nvt}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_em):
        print(f"Error: PDB projection of em file '{args.pdb_proj_em}' does not exist.")
        sys.exit(1)
    
    if not os.path.exists(args.af_temperature_traj):
        print(f"Error: AF temperature of trajectory file '{args.af_temperature_traj}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.af_proj_traj):
        print(f"Error: AF projection of trajectory file '{args.af_proj_traj}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.af_proj_nvt):
        print(f"Error: AF projection of trajectory file '{args.af_proj_nvt}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.af_proj_em):
        print(f"Error: AF projection of trajectory file '{args.af_proj_em}' does not exist.")
        sys.exit(1)

    # Get data
    try:
        print("Reading PDB file...")
        pdb_temperature_time, pdb_temperature_traj = read_xvg(args.pdb_temperature_traj)
        pdb_projection_PC1, pdb_projection_PC2 = read_xvg(args.pdb_proj_traj)
        pdb_nvt_projection_PC1, pdb_nvt_projection_PC2 = read_xvg(args.pdb_proj_nvt)
        pdb_em_projection_PC1, pdb_em_projection_PC2 = read_xvg(args.pdb_proj_em)
        print("Reading AF file...")
        af_temperature_time, af_temperature_traj = read_xvg(args.af_temperature_traj)
        af_projection_PC1, af_projection_PC2 = read_xvg(args.af_proj_traj)
        af_nvt_projection_PC1, af_nvt_projection_PC2 = read_xvg(args.af_proj_nvt)
        af_em_projection_PC1, af_em_projection_PC2 = read_xvg(args.af_proj_em)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)
    
    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)

    #Calculate average temperature for PDB and AF trajectories
    print("\nCalculating average temperatures...")
    all_temperature = np.concatenate([pdb_temperature_traj, af_temperature_traj])
    temp_mean = np.mean(all_temperature)
    print(f"All temperature average: {temp_mean:.2f}")

    #Define constants k_B
    # Boltzmann constant in kJ/(mol·K)
    print("Defining constants...")
    k_B = 0.008314  # kJ/(mol·K)
    
    print("\nGenerating plots...")

    #find max and min in all x(pc1) and y(pc2)
    all_PC1 = np.concatenate([pdb_projection_PC1, pdb_nvt_projection_PC1, pdb_em_projection_PC1, af_projection_PC1, af_nvt_projection_PC1, af_em_projection_PC1])
    all_PC2 = np.concatenate([pdb_projection_PC2, pdb_nvt_projection_PC2, pdb_em_projection_PC2, af_projection_PC2, af_nvt_projection_PC2, af_em_projection_PC2])


    #plot figure
    print(f"Type of af_projection_PC1: {type(af_projection_PC1)}")
    print(f"Length of af_projection_PC1: {len(af_projection_PC1) if hasattr(af_projection_PC1, '__len__') else 'N/A'}")
    print(f"af_projection_PC1 data: {af_projection_PC1}")

    print(f"Type of af_projection_PC2: {type(af_projection_PC2)}")
    print(f"Length of af_projection_PC2: {len(af_projection_PC2) if hasattr(af_projection_PC2, '__len__') else 'N/A'}")
    print(f"af_projection_PC2 data: {af_projection_PC2}")

    # Plot free energy surface for PDB and AF trajectories
    all_F, all_X, all_Y = plot_free_energy_surface(
        all_PC1, all_PC2, 
        temp_mean, k_B, "PC1_PC2", args.save_dir,
        pdb_nvt_projection_PC1, pdb_nvt_projection_PC2,
        pdb_em_projection_PC1, pdb_em_projection_PC2,
        af_nvt_projection_PC1, af_nvt_projection_PC2,
        af_em_projection_PC1, af_em_projection_PC2,
        args.save_time,
        True
    )

if __name__ == "__main__":
    main()


