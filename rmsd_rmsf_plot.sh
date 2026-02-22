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
af_rmsd_dir=""${af_analysis_dir}"/RMSD"
af_rmsf_dir=""${af_analysis_dir}"/RMSF"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_rmsd_dir=""${pdb_analysis_dir}"/RMSD"
pdb_rmsf_dir=""${pdb_analysis_dir}"/RMSF"

#define plot_dir
plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
plot_RMSD_dir="${plot_dir}"/RMSD
plot_RMSF_dir="${plot_dir}"/RMSF


#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$plot_RMSD_dir"
create_dir_if_not_exists "$plot_RMSF_dir"


#define xvg files
pdb_rmsd_file=""${pdb_rmsd_dir}"/traj_internal_rmsd.xvg"
af_rmsd_file=""${af_rmsd_dir}"/traj_fitpdb_rmsd.xvg"

pdb_rmsf_file=""${pdb_rmsf_dir}"/traj_internal_rmsf.xvg"
af_rmsf_file=""${af_rmsf_dir}"/traj_fitpdb_rmsf.xvg"

#set environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py36_env
which python
python_path=$(which python)
expected_path="/home/duanqi/miniconda3/envs/py36_env/bin/python"
if [[ "$python_path" != "$expected_path" ]]; then
    echo "Error: python command is not from the specified py36_env environment"
    echo "Current python path: $python_path"
    echo "Expect python path: $expected_path"
    exit 1
fi

python /home/duanqi/dq/scripts/rmsd_plot.py "${pdb_rmsd_file}" "${af_rmsd_file}" "${plot_RMSD_dir}" "${save_time}"
python /home/duanqi/dq/scripts/rmsf_plot.py "${pdb_rmsf_file}" "${af_rmsf_file}" "${plot_RMSF_dir}" "${save_time}"