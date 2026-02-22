#!/home/duanqi/miniconda3/envs/py36_env/bin/python
# -*- coding: utf-8 -*-
"""
2D轨迹投影对比脚本 (2D Trajectory Projection Comparison Script)

本脚本用于对比PDB和AlphaFold (AF) 结构的主成分分析(PCA)投影结果。
它读取两套投影数据文件(PDB和AF各4个)，在同一个2D图上绘制轨迹散点和参考结构点，
用于直观比较两种方法的采样空间和结构预测的相似性。

功能说明：
  1. 读取XVG格式的投影文件（PC1和PC2坐标）
  2. 处理PDB轨迹、NVT平衡、能量最小化、初始结构的投影
  3. 处理AlphaFold预测结构的相同四类投影
  4. 在统一的坐标系中绘制所有数据
  5. 自动计算坐标轴范围和刻度间隔

输入文件（XVG格式，包含PC1和PC2两列数据）：
  ===== PDB投影文件 (4个) =====
  - pdb_proj_traj    : PDB MD轨迹在PDB-PCA上的投影
  - pdb_proj_nvt     : PDB NVT平衡阶段的投影
  - pdb_proj_em      : PDB能量最小化结构的投影
  - pdb_proj_init    : PDB初始结构的投影
  
  ===== AlphaFold投影文件 (4个) =====
  - af_proj_traj     : AlphaFold预测轨迹（若有）在PDB-PCA上的投影
  - af_proj_nvt      : AlphaFold NVT数据投影
  - af_proj_em       : AlphaFold能量最小化投影
  - af_proj_init     : AlphaFold初始结构投影

输出文件：
  - {save_dir}/proj_2traj_plot_{save_time}.png
    包含：PDB轨迹、AF轨迹、初始结构点的2D散点图

绘图参数：
  - 图片尺寸     : 10" x 6"
  - 分辨率       : 300 DPI
  - 标题         : 'Projection on PDB Principal Components'
  - X轴标签      : 'Projection on PC1(nm)'
  - Y轴标签      : 'Projection on PC2(nm)'
  - 标题字体     : 16pt
  - 轴标签字体   : 14pt
  - 图例字体     : 12pt
  
  数据绘制样式：
    • PDB轨迹      : 黑色小点 (size=2.0, alpha=0.7)
    • AF轨迹       : 蓝色小点 (size=2.0, alpha=0.7)
    • PDB初始结构  : 红色大圆 (size=10.0, 黑色边框)
    • AF初始结构   : 橙色大圆 (size=10.0, 黑色边框)
  
  坐标轴：
    • 自动计算数据范围
    • X/Y轴最多显示10个刻度
    • 刻度间隔向上取整以保证整数

命令行参数：
  python projection_2D2traj.py \
    <pdb_proj_traj> <pdb_proj_nvt> <pdb_proj_em> <pdb_proj_init> \
    <af_proj_traj> <af_proj_nvt> <af_proj_em> <af_proj_init> \
    <save_dir> <save_time>
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
    parser.add_argument('pdb_proj_traj', type=str,
                        help='Path to PDB projection of trajectory file')
    
    parser.add_argument('pdb_proj_nvt', type=str,
                    help='Path to PDB projection of nvt file')
    
    parser.add_argument('pdb_proj_em', type=str,
                help='Path to PDB projection of em file')
    
    parser.add_argument('pdb_proj_init', type=str,
                help='Path to PDB projection of initial structure file')
    
    parser.add_argument('af_proj_traj', type=str,
                        help='Path to AF projection of trajectory file')

    parser.add_argument('af_proj_nvt', type=str,
                        help='Path to AF projection of nvt file')

    parser.add_argument('af_proj_em', type=str,
                help='Path to AF projection of em file') 

    parser.add_argument('af_proj_init', type=str,   
                help='Path to AF projection of initial structure file')  
    
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
    print(f"  PDB projection of trajectory file: {args.pdb_proj_traj}")
    print(f"  PDB projection of nvt file: {args.pdb_proj_nvt}")
    print(f"  PDB projection of em file: {args.pdb_proj_em}")
    print(f"  PDB projection of initial structure file: {args.pdb_proj_init}")
    print(f"  AF projection of trajectory file: {args.af_proj_traj}")
    print(f"  AF projection of nvt file: {args.af_proj_nvt}")
    print(f"  AF projection of em file: {args.af_proj_em}")
    print(f"  AF projection of initial structure file: {args.af_proj_init}")
    print(f"  Save directory: {args.save_dir}")
    
    # Check if the files exist
    if not os.path.exists(args.pdb_proj_traj):
        print(f"Error: PDB projection of trajectory file '{args.pdb_proj_traj}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_nvt):
        print(f"Error: PDB projection of nvt file '{args.pdb_proj_nvt}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_em):
        print(f"Error: PDB projection of em file '{args.pdb_proj_em}' does not exist.")
        sys.exit(1)

    if not os.path.exists(args.pdb_proj_init):
        print(f"Error: PDB projection of initial structure file '{args.pdb_proj_init}' does not exist.")
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

    if not os.path.exists(args.af_proj_init):
        print(f"Error: AF projection of initial structure file '{args.af_proj_init}' does not exist.")
        sys.exit(1)

    # Get data
    try:
        print("Reading PDB file...")
        pdb_projection_PC1, pdb_projection_PC2 = read_xvg(args.pdb_proj_traj)
        pdb_nvt_projection_PC1, pdb_nvt_projection_PC2 = read_xvg(args.pdb_proj_nvt)
        pdb_em_projection_PC1, pdb_em_projection_PC2 = read_xvg(args.pdb_proj_em)
        pdb_init_projection_PC1, pdb_init_projection_PC2 = read_xvg(args.pdb_proj_init)
        print("Reading AF file...")
        af_projection_PC1, af_projection_PC2 = read_xvg(args.af_proj_traj)
        af_nvt_projection_PC1, af_nvt_projection_PC2 = read_xvg(args.af_proj_nvt)
        af_em_projection_PC1, af_em_projection_PC2 = read_xvg(args.af_proj_em)
        af_init_projection_PC1, af_init_projection_PC2 = read_xvg(args.af_proj_init)
    except Exception as e:
        print(f"Error reading XVG files: {str(e)}")
        sys.exit(1)
    
    #Create save directory if it doesn't exist
    os.makedirs(args.save_dir, exist_ok=True)
    
    print("\nGenerating plots...")

    #find max and min in all x(pc1) and y(pc2)
    all_PC1 = np.concatenate([pdb_projection_PC1, pdb_nvt_projection_PC1, pdb_em_projection_PC1, pdb_init_projection_PC1, af_projection_PC1, af_nvt_projection_PC1, af_em_projection_PC1, af_init_projection_PC1])
    all_PC2 = np.concatenate([pdb_projection_PC2, pdb_nvt_projection_PC2, pdb_em_projection_PC2, pdb_init_projection_PC2, af_projection_PC2, af_nvt_projection_PC2, af_em_projection_PC2, af_init_projection_PC2])

    #Use downward rounding to get integer minima and upward rounding to get integer maxima
    x_min = np.floor(np.min(all_PC1))
    x_max = np.ceil(np.max(all_PC1))
    y_min = np.floor(np.min(all_PC2))
    y_max = np.ceil(np.max(all_PC2))

    #plot figure
    print(f"Type of af_projection_PC1: {type(af_projection_PC1)}")
    print(f"Length of af_projection_PC1: {len(af_projection_PC1) if hasattr(af_projection_PC1, '__len__') else 'N/A'}")
    print(f"af_projection_PC1 data: {af_projection_PC1}")

    print(f"Type of af_projection_PC2: {type(af_projection_PC2)}")
    print(f"Length of af_projection_PC2: {len(af_projection_PC2) if hasattr(af_projection_PC2, '__len__') else 'N/A'}")
    print(f"af_projection_PC2 data: {af_projection_PC2}")

    plt.figure(figsize=(10, 6))

    #plot pdb data
    plt.plot(pdb_projection_PC1, pdb_projection_PC2, '.', color='black', markersize=2.0, label='PDB trajectory Projection', alpha=0.7)
    #plot af data
    plt.plot(af_projection_PC1, af_projection_PC2, '.', color='blue', markersize=2.0, label='AlphaFold Projection', alpha=0.7)

    #plot pdb init structure
    plt.plot(pdb_init_projection_PC1, pdb_init_projection_PC2, 'o', color='red', markersize=10.0, alpha=1.0, markeredgecolor='black', markeredgewidth=0.5, label='PDB initial structure')
    #plot af init structure
    plt.plot(af_init_projection_PC1, af_init_projection_PC2, 'o', color='orange', markersize=10.0, alpha=1.0, markeredgecolor='black', markeredgewidth=0.5, label='AF initial structure')





    
    plt.title('Projection on PDB Principal Components', fontsize=16)
    plt.xlabel('Projection on PC1(nm)', fontsize=14)
    plt.ylabel('Projection on PC2(nm)', fontsize=14)

    #Setting the axis range

    #plot legend
    plt.legend(fontsize=12)

    #set x axis ticks
    num_ticks = 10
    step=max(1, int(np.round((x_max - x_min) / (num_ticks - 1))))
    num_steps = (x_max - x_min) // step
    x_max_adjusted = x_min + num_steps * step
    if x_max_adjusted < x_max:
        x_max_adjusted += step

    xticks = np.arange(x_min, x_max_adjusted + step, step)
    plt.xticks(xticks)

    #Setting x axis range
    plt.xlim(x_min, x_max_adjusted)

    #set y axis ticks
    num_ticks = 10
    step=max(1, int(np.round((y_max - y_min) / (num_ticks - 1))))
    num_steps = (y_max - y_min) // step
    y_max_adjusted = y_min + num_steps * step
    if y_max_adjusted < y_max:
        y_max_adjusted += step

    yticks = np.arange(y_min, y_max_adjusted + step, step)
    plt.yticks(yticks)

    #Setting y axis range
    plt.ylim(y_min, y_max_adjusted)

    plt.tight_layout()
    
    #Setting plot path
    plot_path = os.path.join(args.save_dir, f'proj_2traj_plot_{args.save_time}.png')
    plt.savefig(plot_path, dpi=300)
    print(f"Saved 2traj-projection plot to {plot_path}")
if __name__ == "__main__":
    main()


