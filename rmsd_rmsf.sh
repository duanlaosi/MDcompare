#!/bin/bash

#this script aims to 1 calculate rmsd and rmsf 

#set environment
export PATH=/home/duanqi/miniconda3/envs/Gromacs2019.5/bin:$PATH
source ~/miniconda3/etc/profile.d/conda.sh
conda activate Gromacs2019.5
which gmx
gmx_path=$(which gmx)
expected_path="/home/duanqi/miniconda3/envs/Gromacs2019.5/bin/gmx"
if [[ "$gmx_path" != "$expected_path" ]]; then
    echo "Error: gmx command is not from the specified Gromacs2019.5 environment"
    echo "Current gmx path: $gmx_path"
    echo "Expect gmx path: $expected_path"
    export PATH="/home/duanqi/miniconda3/envs/Gromacs2019.5/bin:$PATH"
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
af_rmsd_dir=""${af_analysis_dir}"/RMSD"
af_rmsf_dir=""${af_analysis_dir}"/RMSF"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_PCA_dir=""${pdb_analysis_dir}"/PCA"
pdb_rmsd_dir=""${pdb_analysis_dir}"/RMSD"
pdb_rmsf_dir=""${pdb_analysis_dir}"/RMSF"

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


#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$af_rmsd_dir"
create_dir_if_not_exists "$pdb_rmsd_dir"
create_dir_if_not_exists "$af_rmsf_dir"
create_dir_if_not_exists "$pdb_rmsf_dir"


#----------------------------------------------------------fit af MD traj to af em structure--------------------------------------------------------------------#
#prepare trajectory for rmsd analysis fit af traj to af em structure
#preserve protein only
echo 1 | gmx trjconv -f "$af_results_dir_gromacs"/step_merged.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/traj_af_stripped.xtc

#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
#remove PBC effect and centering the protein
#choose protein1) for centering ; protein1) for output, generate traj_center.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_af_stripped.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/traj_af_mol.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_af_mol.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/traj_af_mol_rt.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/traj_af_mol_rt.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/traj_af_mol_rt_cen.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1


#--------------------------------------------------------fit nvt equilibration traj--------------------------------------------------------------#
#prepare equilibration trajectory for rmsd analysis: fit af traj to af em structure, fit pdb traj to pdb em structure

#  ----------------------------------------AF----------------------------------------  #
# 1 fit to af em structure
#preserve protein only
echo 1 | gmx trjconv -f "$af_results_dir_gromacs"/step4.1_equilibration.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/nvt_af_stripped.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_af_stripped.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/nvt_af_mol.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_af_mol.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/nvt_af_mol_rt.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_af_mol_rt.xtc -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/nvt_af_mol_rt_cen.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

# 2 fit to pdb em structure
#preserve protein only
echo 1 | gmx trjconv -f "$af_results_dir_gromacs"/step4.1_equilibration.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -o "$af_analysis_dir"/nvt_stripped.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_stripped.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/nvt_mol.xtc
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_mol.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$af_analysis_dir"/nvt_mol_rt.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1
echo 1 1 |gmx trjconv -f "$af_analysis_dir"/nvt_mol_rt.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$af_analysis_dir"/nvt_mol_rt_cen.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1

#------------------------------------------PDB-----------------------------------------#
# 1 fit to pdb em structure
#preserve protein only
echo 1 | gmx trjconv -f "$pdb_results_dir_gromacs"/step4.1_equilibration.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -o "$pdb_analysis_dir"/nvt_stripped.xtc
echo 1 1 |gmx trjconv -f "$pdb_analysis_dir"/nvt_stripped.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$pdb_analysis_dir"/nvt_mol.xtc
echo 1 1 |gmx trjconv -f "$pdb_analysis_dir"/nvt_mol.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -fit rot+trans -o "$pdb_analysis_dir"/nvt_mol_rt.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1
echo 1 1 |gmx trjconv -f "$pdb_analysis_dir"/nvt_mol_rt.xtc -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -pbc mol -center yes -o "$pdb_analysis_dir"/nvt_mol_rt_cen.xtc -n "$pdb_analysis_dir"/analysis_index.ndx -tu ns -dt 0.1
#---------------------------------------------------------------------------------------------------------------------#


#-----------------------------------------------------------RMSD using alpha C-------------------------------------------------------------------------------#


# ------------------------MD traj--------------------------- #
#calculate af traj rmsd(use af em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/traj_af_mol_rt_cen.xtc -o "$af_rmsd_dir"/traj_internal_rmsd.xvg -tu ns
#calculate af traj rmsd(use pdb em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/traj_mol_rt_cen.xtc -o "$af_rmsd_dir"/traj_fitpdb_rmsd.xvg -tu ns


#calculate pdb traj rmsd(use pdb em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$pdb_analysis_dir"/traj_mol_rt_cen.xtc -o "$pdb_rmsd_dir"/traj_internal_rmsd.xvg -tu ns


# ------------------nvt equilibration traj------------------ #
#calculate af nvt rmsd(use af em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/nvt_af_mol_rt_cen.xtc -o "$af_rmsd_dir"/nvt_internal_rmsd.xvg -tu ns
#calculate af nvt rmsd(use pdb em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/nvt_mol_rt_cen.xtc -o "$af_rmsd_dir"/nvt_fitpdb_rmsd.xvg -tu ns


#calculate pdb nvt rmsd(use pdb em structure as inderence)
echo C-alpha C-alpha |gmx rms -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$pdb_analysis_dir"/nvt_mol_rt_cen.xtc -o "$pdb_rmsd_dir"/nvt_internal_rmsd.xvg -tu ns


#------------------------------------------------------RMSF using alpha C-----------------------------------------------------#



#calculate af traj rmsf(use af em structure as inderence)
echo C-alpha |gmx rmsf -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/traj_af_mol_rt_cen.xtc -o "$af_rmsf_dir"/traj_internal_rmsf.xvg -res
#calculate af traj rmsf(use pdb em structure as inderence)
echo C-alpha |gmx rmsf -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/traj_mol_rt_cen.xtc -o "$af_rmsf_dir"/traj_fitpdb_rmsf.xvg -res

#calculate pdb traj rmsf(use pdb em structure as inderence)
echo C-alpha |gmx rmsf -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$pdb_analysis_dir"/traj_mol_rt_cen.xtc -o "$pdb_rmsf_dir"/traj_internal_rmsf.xvg -res


#calculate af nvt rmsf(use af em structure as inderence)
echo C-alpha |gmx rmsf -s "$af_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/nvt_af_mol_rt_cen.xtc -o "$af_rmsf_dir"/nvt_internal_rmsf.xvg -res
#calculate af nvt rmsf(use pdb em structure as inderence)
echo C-alpha |gmx rmsf -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$af_analysis_dir"/nvt_mol_rt_cen.xtc -o "$af_rmsf_dir"/nvt_fitpdb_rmsf.xvg -res


#calculate pdb nvt rmsf(use pdb em structure as inderence)
echo C-alpha |gmx rmsf -s "$pdb_results_dir_gromacs"/step4.1_equilibration.tpr -f "$pdb_analysis_dir"/nvt_mol_rt_cen.xtc -o "$pdb_rmsf_dir"/nvt_internal_rmsf.xvg -res
#-----------------------------------------------------------------------------------------------------------#
