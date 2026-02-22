#!/bin/bash

#this script aims to calucalate energy parameter

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
af_energy_dir=""${af_analysis_dir}"/energy"


#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_energy_dir=""${pdb_analysis_dir}"/energy"



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
create_dir_if_not_exists "$af_energy_dir"
create_dir_if_not_exists "$pdb_energy_dir"

#Temperature analysis of pdb and af in MD traj
echo Temperature | gmx energy -f "${pdb_results_dir_gromacs}"/step_combined.edr -o "${pdb_energy_dir}"/traj_temperature.xvg
echo Temperature | gmx energy -f "${af_results_dir_gromacs}"/step_combined.edr -o "${af_energy_dir}"/traj_temperature.xvg

#Pressure analysis of pdb and af in MD traj
echo Pressure | gmx energy -f "${pdb_results_dir_gromacs}"/step_combined.edr -o "${pdb_energy_dir}"/traj_pressure.xvg
echo Pressure | gmx energy -f "${af_results_dir_gromacs}"/step_combined.edr -o "${af_energy_dir}"/traj_pressure.xvg

#Volume analysis of pdb and af in MD traj
echo Volume | gmx energy -f "${pdb_results_dir_gromacs}"/step_combined.edr -o "${pdb_energy_dir}"/traj_volume.xvg
echo Volume | gmx energy -f "${af_results_dir_gromacs}"/step_combined.edr -o "${af_energy_dir}"/traj_volume.xvg

#Total-Energy
echo Total-Energy | gmx energy -f "${pdb_results_dir_gromacs}"/step_combined.edr -o "${pdb_energy_dir}"/traj_total-energy.xvg
echo Total-Energy | gmx energy -f "${af_results_dir_gromacs}"/step_combined.edr -o "${af_energy_dir}"/traj_total-energy.xvg
