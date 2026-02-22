#!/home/duanqi/miniconda3/envs/py36_env/bin/python
# -*- coding: utf-8 -*-
"""
累计方差分析脚本（Cumulative Variance Analysis Script）

本脚本用于分析主成分分析（PCA）的结果。它从XVG文件读取特征值数据，
计算累计方差百分比，并根据指定的阈值百分比选择所需的主成分数量。
最后，生成一张图像展示累计方差曲线、个别方差曲线和阈值线。

功能描述：
  1. read_xvg(filename)          ：从XVG文件读取数据，跳过注释行。
  2. parse_arguments()           ：解析命令行参数。
  3. calculate_cumulative_pc_variance() ：计算累计方差，选择达到阈值的PC。
  4. plot_variance_threshold()   ：绘制累计方差图、个别方差图和阈值线。
  5. main()                      ：主函数，协调整个流程。

命令行参数：
  pdb_eigenvalue         ：PCA特征值文件路径（XVG格式）
  Percentage             ：累计方差阈值百分比（整数，范围0-100）
  save_dir               ：输出图像保存目录

输入文件：
  - XVG格式的特征值文件，通常由GROMACS的gmx covar命令生成。

输出文件：
  - Cumulative_variance.png ：可视化图像，包含：
      * 蓝色曲线 ：累计方差百分比（主要指标）
      * 绿色虚线 ：各PC的个别方差百分比
      * 红色虚线 ：用户指定的阈值百分比
      * 标注点   ：显示具体百分比数值
"""

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
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
    parser.add_argument('pdb_eigenvalue', type=str,
                        help='Path to PDB projection of eigenvalue file')

    parser.add_argument('Percentage', type=int,
                help='PC variance Percentage')   
    
    parser.add_argument('save_dir', type=str,
                        help='Directory to save the plots')
    
    return parser.parse_args()


def calculate_cumulative_pc_variance(eigenvalue_numbers, PC_variances, threshold_percentage):
    """
    Calculate the cumulative PC variance until a specified threshold_percentage is reached
    
    Arguments:
    eigenvalue_numbers: eigenvalue number list
    PC_variances: PC variance list
    threshold_percentage: target pencentage threshold
    
    Return:
    eigenvalue_numbers_filtered: filtered eigenvalue number
    PC_variances_filtered: filtered PC variance
    """
    
    # Calculate total variance summary
    total_variance = np.sum(PC_variances)
    
    # Initialize cumulative variance and index
    cumulative_variance = 0.0
    selected_indices = []
    
    # Traverse each PC variance and calculate the cumulative variance
    for i, pc_variance in enumerate(PC_variances):
        cumulative_variance += pc_variance
        selected_indices.append(i)
        
        # Calculate the cumulative variance as a percentage of the total variance
        cumulative_percentage = (cumulative_variance / total_variance) * 100
        
        # If the target threshold_percentage is reached or exceeded, stop
        if cumulative_percentage >= threshold_percentage:
            break
    
    # Use cycle Extract the corresponding eigenvalue number and PC variance
    eigenvalue_numbers_filtered = [eigenvalue_numbers[i] for i in selected_indices]
    PC_variances_filtered = [PC_variances[i] for i in selected_indices]
    
    return eigenvalue_numbers_filtered, PC_variances_filtered

def plot_variance_threshold(PC_numbers, PC_plot_variances, PC_variances, threshold_percentage, save_dir):
    """
    Plot the cumulative variance and mark the threshold percentage
    
    Arguments:
    PC_numbers: PC number list
    PC_variances: PC variance list
    threshold_percentage: target pencentage threshold
    """
    
    # Calculate cumulative variance
    cumulative_variance = np.cumsum(PC_plot_variances)
    
    # Calculate total variance
    total_variance = np.sum(PC_variances)
    
    # Calculate cumulative percentage of the total variance
    cumulative_percentage = (cumulative_variance / total_variance) * 100

    #Culculate individual percentage
    individual_percentage = (PC_plot_variances / total_variance) * 100

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot cumulative variance
    plt.plot(PC_numbers, cumulative_percentage, marker='o', label='Cumulative Variance')

    # Plot individual variances
    plt.plot(PC_numbers, individual_percentage, marker='s', linestyle='--', label='Individual Variance', color='green')

    # Add a horizontal line for the threshold percentage
    plt.axhline(y=threshold_percentage, color='r', linestyle='--', label=f'Threshold ({threshold_percentage}%)')
    
    # Add labels and title
    plt.xlabel('PC Number')
    plt.ylabel('Percentage (%)')
    plt.title('Cumulative Percentage and individual Percentage of PC')

    # Add xticks
    plt.xticks(PC_numbers)

    # Add grid and legengd
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    # Annotate cumulative percentage points
    for i, val in enumerate(cumulative_percentage):
        plt.annotate(f"{val:.2f}%", 
                    (PC_numbers[i], val), 
                    textcoords="offset points", xytext=(0,8), ha='center', fontsize=8, color='blue')

    # Annotate individual percentage points
    for i, val in enumerate(individual_percentage):
        plt.annotate(f"{val:.2f}%", 
                    (PC_numbers[i], val), 
                    textcoords="offset points", xytext=(0,-12), ha='center', fontsize=8, color='green')

    # Save the plot
    plot_path = os.path.join(save_dir, 'Cumulative_variance.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Saved 2traj-projection plot to {plot_path}")


def main():
    # Parse command line arguments
    args = parse_arguments()
    # Print parameters
    print(f"Running analysis with the following parameters:")
    print(f"  PDB eigenvalue file: {args.pdb_eigenvalue}")
    print(f"  PDB eigenvalue percentage: {args.Percentage}")
    print(f"  Save directory(pdb_PCA_Cumulative-var_dir): {args.save_dir}")
    # Get data
    try:
        print("Reading PDB eigenvalue file...")
        PC_numbers, PC_variances = read_xvg(args.pdb_eigenvalue)

    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)
    
    # Select PCs (cumulatie variance reachng threshold)
    PC_percentage_numbers, PC_percentage_variances = calculate_cumulative_pc_variance(PC_numbers, PC_variances, args.Percentage)
    
    # Calculate PC number (Selected PCs)
    PC_plot_count = len(PC_percentage_numbers)

    # Print result information
    print(f"Choosing {len(PC_percentage_numbers)} PCs")
    print(f"Cumulative Variance Percentage: {(np.sum(PC_percentage_variances) / np.sum(PC_variances)) * 100:.2f}%")
    print(f"PC_{args.Percentage}_numbers: {PC_percentage_numbers}")
    print(f"PC_{args.Percentage}_variances: {PC_percentage_variances}")

    #plot PC_plot_count+1 numbers of PCs (selected PCs and 3+1 behind PCs)
    PC_plot_numbers = PC_numbers[:PC_plot_count + 3]
    PC_plot_variances = PC_variances[:PC_plot_count + 3]
    plot_variance_threshold(PC_plot_numbers, PC_plot_variances, PC_variances, args.Percentage, args.save_dir)
    
if __name__ == "__main__":
    main()