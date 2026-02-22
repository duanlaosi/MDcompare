#!/bin/bash

#this script aims to make movie for pdb trajectory and af trajectory




#define parameter
protein=${1}
af_structure=${2}
af_md_time=${3}

pdb_structure=${4}
pdb_md_time=${5}

save_time=${6}
pdb_colorID=${7}
af_colorID=${8}

#define af directory
af_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${af_structure}"_"${af_md_time}"/"${protein}"_"${af_structure}"_md"
af_results_dir=""${af_basic_dir}"/results"
af_results_dir_gromacs=""${af_results_dir}"/gromacs"
af_analysis_dir=""${af_basic_dir}"/analysis"
af_movie_dir=""${af_analysis_dir}"/movie"

#define pdb directory
pdb_basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${pdb_structure}"_"${pdb_md_time}"/"${protein}"_"${pdb_structure}"_md"
pdb_results_dir=""${pdb_basic_dir}"/results"
pdb_results_dir_gromacs=""${pdb_results_dir}"/gromacs"
pdb_analysis_dir=""${pdb_basic_dir}"/analysis"
pdb_movie_dir=""${pdb_analysis_dir}"/movie"

#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}
create_dir_if_not_exists "$pdb_movie_dir"
create_dir_if_not_exists "$af_movie_dir"

pdb_traj_file="$pdb_analysis_dir"/traj_mol_rt_cen.xtc
pdb_gro_file="$pdb_analysis_dir"/frame0.gro
pdb_colorID=5
af_traj_file="$af_analysis_dir"/traj_mol_rt_cen.xtc
af_gro_file="$af_analysis_dir"/frame0.gro
af_colorID=0

#define f()
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: file $1 not exist!"
        exit 1
    fi
}
check_file "$pdb_traj_file"
check_file "$pdb_gro_file"
check_file "$af_traj_file"
check_file "$af_gro_file"



#make movie for pdb trajectory
echo "Making movie for pdb trajectory..."
bash /home/duanqi/dq/scripts/vmd_movie.tcl "$pdb_gro_file" "$pdb_traj_file" "$save_time" "$pdb_colorID" "$pdb_movie_dir"


bash /home/duanqi/dq/scripts/vmd_movie.tcl "$af_gro_file" "$af_traj_file" "$save_time" "$af_colorID" "$af_movie_dir"

