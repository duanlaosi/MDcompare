#!/bin/bash

#this script aims to make image for pdb and af structure
#set directory and parameters




#define parameter
protein=${1}
af_structure=${2}
af_md_time=${3}
af_color_ID=${4}

pdb_structure=${5}
pdb_md_time=${6}
pdb_color_ID=${7}

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



#define af directory
af_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${af_structure}"_"${af_md_time}"/"${protein}"_"${af_structure}"_md"
af_results_dir=""${af_basic_dir}"/results"
af_results_dir_gromacs=""${af_results_dir}"/gromacs"
af_analysis_dir=""${af_basic_dir}"/analysis"
af_structure_dir=""${af_analysis_dir}"/structure"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_structure_dir=""${pdb_analysis_dir}"/structure"

#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$pdb_structure_dir"
create_dir_if_not_exists "$af_structure_dir"


pdb_gro_file="$pdb_results_dir_gromacs"/input_addh.gro
af_gro_file="$af_results_dir_gromacs"/input_align.gro


#check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: file $1 not exist!"
        exit 1
    fi
}
check_file "$pdb_gro_file"
check_file "$af_gro_file"



#make image for pdb structure
echo "Making movie for pdb trajectory..."
bash /home/duanqi/dq/scripts/vmd_structure.tcl "$pdb_gro_file" "$pdb_structure_dir" "$pdb_color_ID"

#make image for af structure
echo "Making movie for af trajectory..."
bash /home/duanqi/dq/scripts/vmd_structure.tcl "$af_gro_file" "$af_structure_dir" "$af_color_ID"

