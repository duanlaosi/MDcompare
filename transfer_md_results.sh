#!/bin/bash

protein=$1
structure=$2
md_time=$3


basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${structure}"_"${md_time}"/"${protein}"_"${structure}"_md"
results_dir=""${basic_dir}"/results"

create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}

create_dir_if_not_exists "$basic_dir"
create_dir_if_not_exists "$results_dir"

#transport gromacs directory to results directory
#get results/gromacs directory
target_dir="/home/duanqi/"${protein}"/"${protein}"_"${structure}"/gromacs"
#copy  simulation results from server 172.31.152.81 to server 172.31.151.118
scp -r duanqi@172.31.152.81:"${target_dir}" duanqi@172.31.151.118:"${results_dir}"/

#sign in 172.31.151.118 to run other scripts
if [ $? -eq 0 ]; then
    echo "successfully transfer"
else
    echo "failed trasfer, end"
    exit 1
fi