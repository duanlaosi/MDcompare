#!/bin/bash
# -*- coding: utf-8 -*-
export LANG=zh_CN.UTF-8
export LC_ALL=zh_CN.UTF-8

# 本脚本用于对PDB轨迹进行主成分分析（PCA），
# 并将轨迹及相关结构投影到主要主成分上。
# 同时对轨迹和结构进行拟合到能量最小化的参考结构，结果保存于分析目录内。
#
# 使用的输入文件包括:
#   - step4.1_equilibration.tpr : 用于生成索引和参考对齐
#   - step_merged.xtc          : MD轨迹文件，随后被处理和投影
#   - 各种结构文件如 step4.0_minimization.gro、
#     step4.1_equilibration.gro、step3_input.gro 等，用于投影
# 所有生成的PCA输出（eigenvalues, eigenvectors, proj_xvg 等）
# 会被移动到分析目录下的PCA子目录（${analysis_dir}/PCA）。
# 可以在原始分析目录中找到中间文件如 traj_mol_rt_cen.xtc 等。

#this script aims to 1 calculate PCs in pdb trajectory and then project on pdb PCs
#this script aims to 2 fit pdb trajectory and structures to pdb minimization structure save in analysis_dir

#set environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate Gromacs2019.5
which gmx
gmx_path=$(which gmx)
expected_path="/home/duanqi/miniconda3/envs/Gromacs2019.5/bin/gmx"
if [[ "$gmx_path" != "$expected_path" ]]; then
    echo "Error: gmx command is not from the specified Gromacs2019.5 environment"
    echo "Current gmx path: $gmx_path"
    echo "Expect gmx path: $expected_path"
    exit 1
fi

#set parameter
protein=$1
structure=$2
md_time=$3

#define pdb directory 
basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${structure}"_"${md_time}"/"${protein}"_"${structure}"_md"
results_dir=""${basic_dir}"/results"
#define results directory and analysis directory
results_dir_gromacs=""${results_dir}"/gromacs"
analysis_dir=""${basic_dir}"/analysis"
PCA_dir=""${analysis_dir}"/PCA"


#make analysis_index from md simulation structure pdbter energy minimization
gmx make_ndx -f "$results_dir_gromacs"/step4.1_equilibration.tpr -o "$analysis_dir"/analysis_index.ndx <<EOF
0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
q
EOF

#prepare trajectory for PCA and other analysis

#preserve protein only
echo 1 | gmx trjconv -f "$results_dir_gromacs"/step_merged.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -o "$analysis_dir"/traj_stripped.xtc
#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
echo 1 1 |gmx trjconv -f "$analysis_dir"/traj_stripped.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$analysis_dir"/traj_mol.xtc

echo 1 1 |gmx trjconv -f "$analysis_dir"/traj_mol.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$analysis_dir"/traj_mol_rt.xtc -n "$analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

echo 1 1 |gmx trjconv -f "$analysis_dir"/traj_mol_rt.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$analysis_dir"/traj_mol_rt_cen.xtc -n "$analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

#calculate the covariance matrix(PCA)

#select "C-alpha"(3) for least-squares fit, select "C-alpha"(3) for PCA analysis
echo 3 3 |gmx covar -f "$analysis_dir"/traj_mol_rt_cen.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -o "$analysis_dir"/traj_eigenval.xvg -v "$analysis_dir"/traj_eigenvec.trr -tu ns -dt 0.1 -av "$analysis_dir"/traj_average.pdb

#projections

#project pdb md trajectory(traj_mol_rt_cen.xtc) on pdb PC1 PC2 (C-alpha)
echo 3 3 |gmx anaeig -f "$analysis_dir"/traj_mol_rt_cen.xtc -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$analysis_dir"/traj_eigenvec.trr -proj "$analysis_dir"/proj_traj.xvg -2d "$analysis_dir"/2Dproj_traj.xvg -tu ns

#preserve protein only 
echo 1 1 | gmx trjconv -f "$results_dir_gromacs"/step4.0_minimization.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$analysis_dir"/step4.0_minimization_stripped.gro
#project pdb em structure(step4.0_minimization.gro) on pdb PC1 PC2
echo 3 3 |gmx anaeig -f "$analysis_dir"/step4.0_minimization_stripped.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$analysis_dir"/traj_eigenvec.trr -proj "$analysis_dir"/proj_em.xvg -2d "$analysis_dir"/2Dproj_em.xvg -tu ns

#preserve protein only
echo 1 1 | gmx trjconv -f "$results_dir_gromacs"/step4.1_equilibration.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$analysis_dir"/step4.1_equilibration_stripped.gro
#project pdb nvt structure(step4.1_equilibration.gro) on pdb PC1 PC2
echo 3 3 |gmx anaeig -f "$analysis_dir"/step4.1_equilibration_stripped.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$analysis_dir"/traj_eigenvec.trr -proj "$analysis_dir"/proj_nvt.xvg -2d "$analysis_dir"/2Dproj_nvt.xvg -tu ns

#preserve protein only
echo 1 1 | gmx trjconv -f "$results_dir_gromacs"/step3_input.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$analysis_dir"/step3_input_stripped.gro
#project pdb initial structure(step3_input.gro) on pdb PC1 PC2
echo 3 3 |gmx anaeig -f "$analysis_dir"/step3_input_stripped.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$analysis_dir"/traj_eigenvec.trr -proj "$analysis_dir"/proj_init.xvg -2d "$analysis_dir"/2Dproj_init.xvg -tu ns

# #extract first frame in trajectory
# echo 1 1 |gmx trjconv -s "$results_dir_gromacs"/step4.1_equilibration.tpr -f "$analysis_dir"/traj_mol_rt_cen.xtc -o "$analysis_dir"/frame0.gro -b 0 -e 0
# #project frame0 structure(frame0.gro)
# echo 3 3 |gmx anaeig -f "$analysis_dir"/frame0.gro -s "$results_dir_gromacs"/step4.1_equilibration.tpr -n "$analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$analysis_dir"/traj_eigenvec.trr -proj "$analysis_dir"/proj_frame0.xvg -2d "$analysis_dir"/2Dproj_frame0.xvg -tu ns

#save PCA results into /$PCA_dir directory
mv "$analysis_dir"/traj_eigenvec.trr "$analysis_dir"/traj_eigenval.xvg "$analysis_dir"/proj* "$analysis_dir"/2Dproj* "${PCA_dir}"/