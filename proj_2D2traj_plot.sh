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

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_PCA_dir=""${pdb_analysis_dir}"/PCA"

#define plot_dir
plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
plot_PCA_dir="${plot_dir}"/PCA


#define traj projection xvg files
pdb_file=""${pdb_PCA_dir}"/2Dproj_traj.xvg"
af_file=""${af_PCA_dir}"/2Dproj_traj.xvg"
#define projection xvg files 
pdb_proj_traj=""${pdb_PCA_dir}"/2Dproj_traj.xvg"
pdb_proj_nvt=""${pdb_PCA_dir}"/2Dproj_nvt.xvg"
pdb_proj_em=""${pdb_PCA_dir}"/2Dproj_nvt.xvg"
pdb_proj_init=""${pdb_PCA_dir}"/2Dproj_init.xvg"
af_proj_traj=""${af_PCA_dir}"/2Dproj_traj.xvg"
af_proj_nvt=""${af_PCA_dir}"/2Dproj_nvt.xvg"
af_proj_em=""${af_PCA_dir}"/2Dproj_em.xvg"
af_proj_init=""${af_PCA_dir}"/2Dproj_init.xvg"


#set environment
export PATH=/home/duanqi/miniconda3/envs/py38_env/bin:$PATH
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

python /home/duanqi/dq/scripts/projection_2D2traj.py "${pdb_proj_traj}" "${pdb_proj_nvt}" "${pdb_proj_em}" "${pdb_proj_init}" "${af_proj_traj}" "${af_proj_nvt}" "${af_proj_em}" "${af_proj_init}" "${plot_PCA_dir}" "${save_time}"
#python /home/duanqi/dq/scripts/interval_test.py "${pdb_file}" "${af_file}" "${plot_dir}" "${save_time}"
