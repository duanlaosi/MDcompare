#!/home/duanqi/miniconda3/envs/py38_env/bin/python
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import os
import sys
import argparse
def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Analysis Script')
    
    # Define positional arguments (required, no default values)
    parser.add_argument('pdb_structure_gro', type=str,
                        help='Path to PDB structure Gro file')

    parser.add_argument('pdb_secondary_structure_dir', type=str,
                        help='Path to PDB structure Directory')
    
    parser.add_argument('af_structure_gro', type=str,
                    help='Path to AF structure Gro file')
    
    parser.add_argument('af_secondary_structure_dir', type=str,
                        help='Path to AF structure Directory')   
    
    return parser.parse_args()

def simple_gro_dssp(gro_file_path):
    """
    使用MDTraj从gro文件计算DSSP二级结构
    
    Args:
        gro_file_path (str): gro文件路径
        
    Returns:
        str: DSSP二级结构字符串
    """
    try:
        # 加载gro文件
        # 注意：gro文件通常需要配合topology文件，这里假设gro文件包含完整信息
        traj = md.load(gro_file_path)
        
        # 计算DSSP
        dssp = md.compute_dssp(traj, simplified=True)
        
        # 获取第一帧的二级结构（如果有多帧）
        dssp_string = ''.join(dssp[0])
        
        return dssp_string
        
    except Exception as e:
        print(f"Error processing {gro_file_path}: {e}")
        return None

def main():
    # Parse command line arguments
    args = parse_arguments()
    pdb_gro_dssp_string = simple_gro_dssp(args.pdb_structure_gro)
    af_gro_dssp_string = simple_gro_dssp(args.af_structure_gro)

if __name__ == "__main__":
    main()