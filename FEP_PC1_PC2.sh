#!/bin/bash

#define parameter
protein=$1
af_structure=$2
af_md_time=$3

pdb_structure=$4
pdb_md_time=$5

save_time=$6

#define af directory
af_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${af_structure}"_"${af_md_time}"/"${protein}"_"${af_structure}"_md"
af_results_dir=""${af_basic_dir}"/results"
af_results_dir_gromacs=""${af_results_dir}"/gromacs"
af_analysis_dir=""${af_basic_dir}"/analysis"
af_PCA_dir=""${af_analysis_dir}"/PCA"
#define energy directory
af_energy_dir=""${af_analysis_dir}"/energy"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_PCA_dir=""${pdb_analysis_dir}"/PCA"
#define energy directory
pdb_energy_dir=""${pdb_analysis_dir}"/energy"

#define plot_dir
plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
plot_FEP_dir="${plot_dir}"/FEP
plot_FEP_PC1_PC2_dir="${plot_FEP_dir}"/FEP_PC1_PC2

#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$plot_FEP_dir"
create_dir_if_not_exists "$plot_FEP_PC1_PC2_dir"

#define projection xvg files 
pdb_proj_traj=""${pdb_PCA_dir}"/2Dproj_traj.xvg"
pdb_proj_nvt=""${pdb_PCA_dir}"/2Dproj_nvt.xvg"
pdb_proj_em=""${pdb_PCA_dir}"/2Dproj_em.xvg"
af_proj_traj=""${af_PCA_dir}"/2Dproj_traj.xvg"
af_proj_nvt=""${af_PCA_dir}"/2Dproj_nvt.xvg"
af_proj_em=""${af_PCA_dir}"/2Dproj_em.xvg"
#define edr  files 
pdb_temperature_traj=""${pdb_energy_dir}"/traj_temperature.xvg"
af_temperature_traj=""${af_energy_dir}"/traj_temperature.xvg"



#set environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py38_env
which python
python_path=$(which python)
expected_path="/home/duanqi/miniconda3/envs/py38_env/bin/python"
if [[ "$python_path" != "$expected_path" ]]; then
    echo "Error: python command is not from the specified py38_env environment"
    echo "Current python path: $python_path"
    echo "Expect python path: $expected_path"
    exit 1
fi


python /home/duanqi/dq/scripts/FEP_PC1_PC2.py "${pdb_temperature_traj}" "${pdb_proj_traj}" "${pdb_proj_nvt}" "${pdb_proj_em}" "${af_temperature_traj}" "${af_proj_traj}" "${af_proj_nvt}" "${af_proj_em}" "${plot_FEP_PC1_PC2_dir}" "${save_time}"

