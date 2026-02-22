#!/home/duanqi/miniconda3/envs/py36_env/bin/python
import matplotlib.pyplot as plt
import numpy as np
from pymbar import timeseries
from pymbar.timeseries import subsampleCorrelatedData
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
    
    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
   #------------------------------------------------------------------------------------------# 
    data = pdb_projection_PC1
    g = timeseries.statisticalInefficiency(data)
    print(f"Statistical inefficiency g ≈ {g:.2f}")
    suggested_interval = int(np.ceil(g))
    print(f"Recommended sampling interval: {suggested_interval} frames")

    indices = timeseries.subsampleCorrelatedData(data)
    decorrelated_data = data[indices]
    def autocorrelation(x):
        x = x - np.mean(x)
        result = np.correlate(x, x, mode='full')
        result = result[result.size // 2:]
        return result / result[0]

    acf = autocorrelation(data)
    plt.plot(acf[:300])
    plt.ylabel("Autocorrelation")
    plt.title("Autocorrelation of PC1 Projection")
    plt.grid()
    plt.show()
    #------------------------------------------------------------------------------------------------------#
    data = af_projection_PC1
    g = timeseries.statisticalInefficiency(data)
    print(f"Statistical inefficiency g ≈ {g:.2f}")
    suggested_interval = int(np.ceil(g))
    print(f"Recommended sampling interval: {suggested_interval} frames")

    indices = timeseries.subsampleCorrelatedData(data)
    decorrelated_data = data[indices]
    def autocorrelation(x):
        x = x - np.mean(x)
        result = np.correlate(x, x, mode='full')
        result = result[result.size // 2:]
        return result / result[0]

    acf = autocorrelation(data)
    plt.plot(acf[:300])
    plt.ylabel("Autocorrelation")
    plt.title("Autocorrelation of PC1 Projection")
    plt.grid()
    plt.show()



if __name__ == "__main__":
    main()


