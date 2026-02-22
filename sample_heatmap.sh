#!/bin/bash

#define parameter
protein=$1

af_structure=$2
af_md_time=$3

pdb_structure=$4
pdb_md_time=$5

frame_number=$6
group_number=$7

min_interval=$8
max_interval=$9
num_intervals=${10}

save_time=${11}

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

#define plot_dir
plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
plot_PCA_dir="${plot_dir}"/PCA
plot_PCA_heatmap_dir="${plot_dir}"/PCA/heatmap
#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$plot_PCA_dir"
create_dir_if_not_exists "$plot_PCA_heatmap_dir"

#define xvg files
pdb_file=""${pdb_PCA_dir}"/2Dproj_traj.xvg"
af_file=""${af_PCA_dir}"/2Dproj_traj.xvg"

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

python /home/duanqi/dq/scripts/sample_heatmap.py "${pdb_file}" "${af_file}" "${frame_number}" "${group_number}" "${min_interval}" "${max_interval}" "${num_intervals}" "${plot_PCA_heatmap_dir}" "${save_time}"
