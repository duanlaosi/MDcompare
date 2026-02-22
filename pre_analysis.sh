#!/bin/bash

#define parameter
protein=$1
structure=$2
md_time=$3
n_step=$4
#define results directory and analysis directory
basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${structure}"_"${md_time}"/"${protein}"_"${structure}"_md"
results_dir=""${basic_dir}"/results"
results_dir_gromacs=""${results_dir}"/gromacs"
analysis_dir=""${basic_dir}"/analysis"
PCA_dir=""${analysis_dir}"/PCA"
rmsd_dir=""${analysis_dir}"/RMSD"
rmsf_dir=""${analysis_dir}"/RMSF"
plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
plot_rmsd_dir="${plot_dir}"/RMSD
plot_rmsf_dir="${plot_dir}"/RMSF
plot_PCA_dir="${plot_dir}"/PCA
plot_PCA_heatmap_dir="${plot_dir}"/PCA/heatmap
#Function to check and create directories if they don't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        mkdir -p "$1"
    fi
}

create_dir_if_not_exists "$basic_dir"
create_dir_if_not_exists "$results_dir"
create_dir_if_not_exists "$results_dir_gromacs"
create_dir_if_not_exists "$analysis_dir"
create_dir_if_not_exists "$PCA_dir"
create_dir_if_not_exists "$rmsd_dir"
create_dir_if_not_exists "$rmsf_dir"
create_dir_if_not_exists "$plot_dir"
create_dir_if_not_exists "$plot_PCA_dir"
create_dir_if_not_exists "$plot_PCA_heatmap_dir"
create_dir_if_not_exists "$plot_rmsd_dir"
create_dir_if_not_exists "$plot_rmsf_dir"

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
    exit 1
fi

#change into result/gromacs directory
cd "$results_dir_gromacs"
#delete log files
rm *.log
rm *.cpt

#generate input files
#choose 0 to step5_1 and c to step5_*
xtc_files=""
edr_files=""
input_list="0"
for ((i=1; i<=n_step; i++)); do
    xtc_files+="step5_${i}.xtc "
    edr_files+="step5_${i}.edr "
    if [[ $i -gt 1 ]]; then
        input_list+="\nc"
    fi
done

#combine xtc files and fix time to make it continuous
echo -e "$input_list" | gmx trjcat -f $xtc_files -o step_merged.xtc -settime yes

#combine edr files and fix time
echo -e "$input_list" | gmx eneconv -f $edr_files -o step_combined.edr -settime yes

#remove step5_1* ~  step5_50*  
rm step5_{1..9}* step5_{10..50}*


#return analysis directory
cd "${analysis_dir}"