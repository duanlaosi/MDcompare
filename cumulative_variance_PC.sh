#!/bin/bash
# -*- coding: utf-8 -*-
export LANG=zh_CN.UTF-8
export LC_ALL=zh_CN.UTF-8

# 本脚本计算PDB主成分分析结果的累计方差，并根据指定阈值
# 选择所需的主成分数量。它调用位于
# /home/duanqi/dq/scripts/cumulative_variance.py 的Python脚本
# 来生成图像。
#
# 使用的输入文件：
#   - ${pdb_PCA_dir}/traj_eigenval.xvg : 来自PCA的特征值文件。
#
# 生成的输出文件位置：
#   - 图像保存到 ${pdb_PCA_Cumulative_var_dir}/Cumulative_variance.png
#     目录由脚本在运行前创建（如果不存在）。
#
# 其他部分保持不变。

#define parameter
protein=$1

pdb_structure=$2
pdb_md_time=$3


Percentage=$4

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_PCA_dir=""${pdb_analysis_dir}"/PCA"
pdb_PCA_Cumulative_var_dir="${pdb_PCA_dir}"/Cumulative_var
#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "${pdb_PCA_Cumulative_var_dir}"

#define xvg files
pdb_file="${pdb_PCA_dir}"/traj_eigenval.xvg 

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

python /home/duanqi/dq/scripts/cumulative_variance.py "${pdb_file}" "${Percentage}" "${pdb_PCA_Cumulative_var_dir}"
