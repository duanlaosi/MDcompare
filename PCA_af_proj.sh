#!/bin/bash

#this script aims to 1 project af trajectory and structures on pdb PCs
#this script aims to 2 fit af trajectory and structures to pdb minimization structure and save in af_analysis_dir

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
af_structure=$2
af_md_time=$3

pdb_structure=$4
pdb_md_time=$5

#define af directory
af_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${af_structure}"_"${af_md_time}"/"${protein}"_"${af_structure}"_md"
af_results_dir=""${af_basic_dir}"/results"
af_results_dir_gromacs=""${af_results_dir}"/gromacs"
af_analysis_dir=""${af_basic_dir}"/analysis"
af_PCA_dir=""${af_analysis_dir}"/PCA"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_PCA_dir=""${pdb_analysis_dir}"/PCA"

#prepare trajectory for PCA and other analysis
#preserve protein only
echo 1 | gmx trjconv -f "$af_results_dir_gromacs"/step_merged.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/traj_stripped.xtc

#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_stripped.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/traj_mol.xtc

echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_mol.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/traj_mol_rt.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_mol_rt.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/traj_mol_rt_cen.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

#projections

#project af md trajectory(traj_mol_rt_cen.xtc) on pdb PC1 PC2 (C-alpha)
echo 3 3 |gmx anaeig -f "$af_analysis_dir"/traj_mol_rt_cen.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -n "$pdb_analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$pdb_PCA_dir"/traj_eigenvec.trr -proj "$af_analysis_dir"/proj_traj.xvg -2d "$af_analysis_dir"/2Dproj_traj.xvg -tu ns

#preserve protein only 
echo 1 1 | gmx trjconv -f "$af_results_dir_gromacs"/step4.0_minimization.gro -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/step4.0_minimization_stripped.gro
echo 1 1 | gmx trjconv -f "$af_analysis_dir"/step4.0_minimization_stripped.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/step4.0_minimization_fit.gro
#project af em structure(step4.0_minimization.gro) on pdb PC1 PC2
echo 3 3 |gmx anaeig -f "$af_analysis_dir"/step4.0_minimization_fit.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -n "$pdb_analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$pdb_PCA_dir"/traj_eigenvec.trr -proj "$af_analysis_dir"/proj_em.xvg -2d "$af_analysis_dir"/2Dproj_em.xvg -tu ns

#preserve protein only
echo 1 1 | gmx trjconv -f "$af_results_dir_gromacs"/step4.1_equilibration.gro -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/step4.1_equilibration_stripped.gro
echo 1 1 | gmx trjconv -f "$af_analysis_dir"/step4.1_equilibration_stripped.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/step4.1_equilibration_fit.gro
#project af nvt structure(step4.1_equilibration.gro) on pdb PC1 PC2(Calpha)
echo 3 3 |gmx anaeig -f "$af_analysis_dir"/step4.1_equilibration_fit.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -n "$pdb_analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$pdb_PCA_dir"/traj_eigenvec.trr -proj "$af_analysis_dir"/proj_nvt.xvg -2d "$af_analysis_dir"/2Dproj_nvt.xvg -tu ns

#preserve protein only
echo 1 1 | gmx trjconv -f "$af_results_dir_gromacs"/step3_input.gro -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/step3_input_stripped.gro
echo 1 1 | gmx trjconv -f "$af_analysis_dir"/step3_input_stripped.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/step3_input_fit.gro
#project af initial structure(step3_input.gro) on pdb PC1 PC2
echo 3 3 |gmx anaeig -f "$af_analysis_dir"/step3_input_fit.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -n "$pdb_analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$pdb_PCA_dir"/traj_eigenvec.trr -proj "$af_analysis_dir"/proj_init.xvg -2d "$af_analysis_dir"/2Dproj_init.xvg -tu ns


# #extract first frame in trajectory
# echo 1 1 |gmx trjconv -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/traj_mol_rt_cen.xtc -o "$af_analysis_dir"/frame0.gro -b 0 -e 0
# #project frame0 structure(frame0.gro)
# echo 3 3 |gmx anaeig -f "$af_analysis_dir"/frame0.gro -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -n "$pdb_analysis_dir"/analysis_index.ndx -first 1 -last 2 -v "$pdb_PCA_dir"/traj_eigenvec.trr -proj "$af_analysis_dir"/proj_frame0.xvg -2d "$af_analysis_dir"/2Dproj_frame0.xvg -tu ns

mv "$af_analysis_dir"/proj* "$af_analysis_dir"/2Dproj* "${af_PCA_dir}"/