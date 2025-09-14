#!/bin/bash

# transfer md results from 172.31.152.81 to 172.31.151.118 and create basic_dir/results | protein name &structure type and MDtime
# basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${structure}"_"${md_time}"/"${protein}"_"${structure}"_md"
# results_dir=""${basic_dir}"/results"
bash /home/duanqi/dq/scripts/transfer_md_results.sh 1BPI pdb0 4_25
bash /home/duanqi/dq/scripts/transfer_md_results.sh 1BPI pdb0 5_5
bash /home/duanqi/dq/scripts/transfer_md_results.sh 1BPI af 5_5
bash /home/duanqi/dq/scripts/transfer_md_results.sh 1BPI pdb_SS 5_18
bash /home/duanqi/dq/scripts/transfer_md_results.sh 1BPI af_SS 5_18
bash /home/duanqi/dq/scripts/transfer_md_results.sh 6ZD3 pdb_triclinic 6_7
bash /home/duanqi/dq/scripts/transfer_md_results.sh 6ZD3 af_triclinic 6_7



# check and create all directories needed(results/gromacs basic_dir/analysis analysis/ /home/duanqi/dq/"${protein}"/plot) |protein name &structure type and MDtime&step number
# basic_dir="/home/duanqi/dq/"${protein}"/"${protein}"_"${structure}"_"${md_time}"/"${protein}"_"${structure}"_md"
# results_dir=/home/duanqi/dq/1BPI/1BPI_pdb_SS_5_18/1BPI_pdb_SS_md/results
# results_dir_gromacs=/home/duanqi/dq/1BPI/1BPI_pdb_SS_5_18/1BPI_pdb_SS_md/results/gromacs
# analysis_dir=/home/duanqi/dq/1BPI/1BPI_pdb_SS_5_18/1BPI_pdb_SS_md/analysis
# plot_dir="/home/duanqi/dq/"${protein}"/plot/${pdb_structure}_${af_structure}"
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI pdb0 4_25
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI af 4_28 
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI pdb0 5_5 50
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI af 5_5 50
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI pdb_SS 5_18 50
bash /home/duanqi/dq/scripts/pre_analysis.sh 1BPI af_SS 5_18 50
bash /home/duanqi/dq/scripts/pre_analysis.sh 6ZD3 pdb_triclinic 6_7 50
bash /home/duanqi/dq/scripts/pre_analysis.sh 6ZD3 af_triclinic 6_7 50



#-------------------------------------------------fit traj and PCA---------------------------------------------------------#



#manipulate pdb trajectory to fit pdb energy minimization structure
#PCA analysis for pdb structure and use pdb structure caluculate covariance matrix and PCs
#save results in /PCA save trajectory and structure in analysis
bash /home/duanqi/dq/scripts/PCA_pdb.sh 1BPI pdb0 4_25
bash /home/duanqi/dq/scripts/PCA_pdb.sh 1BPI pdb0 5_5
bash /home/duanqi/dq/scripts/PCA_pdb.sh 1BPI pdb_SS 5_18
bash /home/duanqi/dq/scripts/PCA_pdb.sh 6ZD3 pdb_triclinic 6_7

# #############PC cumulative variance(py38_env)####################
bash /home/duanqi/dq/scripts/cumulative_variance_PC.sh 6ZD3 pdb_triclinic 6_7 90
bash /home/duanqi/dq/scripts/cumulative_variance_PC.sh 1BPI pdb_SS 5_18 90



#manipulate af trajectpry to fit pdb energy minimization structure
#af structure project on pdb PCs
#save results in /PCA save trajectory and structure in analysis
# project two traj on PDB PC1 PC2 
bash /home/duanqi/dq/scripts/PCA_af_proj.sh 1BPI af 4_28 pdb0 4_25
bash /home/duanqi/dq/scripts/PCA_af_proj.sh 1BPI af 5_5 pdb0 5_5
bash /home/duanqi/dq/scripts/PCA_af_proj.sh 1BPI af_SS 5_18 pdb_SS 5_18
bash /home/duanqi/dq/scripts/PCA_af_proj.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7


#----------------------------------------------------Projection analysis-----------------------------------------------------------#


###projection plot###----------------------(PC1 PC2)

# plot 2traj 2D projection of PC1 and PC2 plot(py38_env)
bash /home/duanqi/dq/scripts/proj_2D2traj_plot.sh 1BPI af 4_28 pdb0 4_25
bash /home/duanqi/dq/scripts/proj_2D2traj_plot.sh 1BPI af 5_5 pdb0 5_5 500ns
bash /home/duanqi/dq/scripts/proj_2D2traj_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/proj_2D2traj_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns

# -----------------------------------------------------PC1 plot-----------------------------------------------------------#
# plot 2traj 1D projection of PC1(py38_env)--------------------------------------------------------------------------------------------------------------------------------#
bash /home/duanqi/dq/scripts/proj_1D2traj_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns
bash /home/duanqi/dq/scripts/proj_1D2traj_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns

# plot 2traj 1D distribution of PC1(py38_env)
bash /home/duanqi/dq/scripts/distribution_1D2traj_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/distribution_1D2traj_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns

###T-TEST and LEVENE-TEST###(PC1)

# plot ttest and levene test results of different sample intervals min max and number of intervals(Using xvvg file's first column to do t-test and levene-test)
bash /home/duanqi/dq/scripts/sample_P_plot.sh 1BPI af 4_28 pdb0 4_25
bash /home/duanqi/dq/scripts/sample_P_plot.sh 1BPI af 5_5 pdb0 5_5 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 1.0 10 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 1.1 2.0 10 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 7.9 10 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 9.0 10 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 10.0 10 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 20.0 20 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 10.0 30.0 20 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 30.0 50.0 20 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 40.0 60.0 20 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 10.0 11 500ns
bash /home/duanqi/dq/scripts/sample_P_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 0.1 1.5 15 500ns

# plot bootstrap sample p-value in different sample interval  (define min_interval & max-interval & PDB-AF)
bash /home/duanqi/dq/scripts/bootstrap_sample_P_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 0.1 1.5 15 500ns
bash /home/duanqi/dq/scripts/bootstrap_sample_P_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 10.0 11 500ns

# # plot 10 group internal p-value in different sample interval  (define min_interval & max-interval & interval-groups) and group number (10)
# bash /home/duanqi/dq/scripts/sample_P_internal.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 0.1 1.5 15 500ns
# bash /home/duanqi/dq/scripts/sample_P_internal.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 10.0 11 500ns

# plot bootstrap sample p-value in different sample interval  (define min_interval & max-interval & interval-groups) and group number (10)
bash /home/duanqi/dq/scripts/bootstrap_sample_P_internal.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 0.1 1.5 15 500ns
bash /home/duanqi/dq/scripts/bootstrap_sample_P_internal.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 0.1 10.0 11 500ns



# # plot 2traj ttest heatmap plot(PdbID af af_time pdb_structure_name pdb_time md_length frame_number group_number)
# bash /home/duanqi/dq/scripts/heat_map_plot.sh 1BPI af 5_5 pdb0 5_5 5001 10 500ns
# bash /home/duanqi/dq/scripts/heat_map_plot.sh 1BPI af 4_28 pdb0 4_25 1001 10 100ns
# bash /home/duanqi/dq/scripts/heat_map_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 500ns

# plot heat map (t test and levene test) in different sample interval  (define min_interval & max-interval & interval-groups) and group number (10)
bash /home/duanqi/dq/scripts/sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 0.1 1.0 10 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 1.6 1.8 3 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 1.9 2.0 2 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 0.6 2.0 15 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 2.1 3.0 2 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 9.0 10.0 2 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 5 10.0 10.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 10.0 10.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 2 1.0 1.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 2 1.3 1.3 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 40.0 40.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 50.0 50.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 30.0 30.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 45.0 45.0 1 500ns
bash /home/duanqi/dq/scripts/sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 1.0 2.0 2 500ns

###ACF###------------------------------(PC1)

# plot ACF of pdb projection and af projection
bash /home/duanqi/dq/scripts/ACF_test.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns

# -----------------------------------------------------PC1 & PC2 plot-----------------------------------------------------------#
# Two PC ANOVA (Hotelling's T-squared test) plot
bash /home/duanqi/dq/scripts/sample_MANOVA.sh 1BPI af_SS 5_18 pdb_SS 5_18 0.1 1.5 15 500ns
bash /home/duanqi/dq/scripts/sample_MANOVA.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 10.0 11 500ns



#-----------------------------------------RMSD and RMSF------------------------------------------------------#


# calculate rmsd and RMSF for pdb and af
# Af traj use pdb em structure as reference and calculate RMSD saved as traj_fitpdb_rmsd.xvg
# Af traj use af em structure as reference and calculate RMSD saved as traj_internal_rmsd.xvg
# PDB traj use pdb em structure as reference and calculate RMSD saved as traj_internal_rmsd.xvg
# Af nvt traj use pdb em structure as reference and calculate RMSD saved as traj_fitpdb_rmsd.xvg
# Af nvt traj use af em structure as reference and calculate RMSD saved as traj_internal_rmsd.xvg
# PDB nvt traj use pdb em structure as reference and calculate RMSD saved as traj_internal_rmsd.xvg
bash /home/duanqi/dq/scripts/rmsd_rmsf.sh 1BPI af_SS 5_18 pdb_SS 5_18
bash /home/duanqi/dq/scripts/rmsd_rmsf.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7


#---------------------------------------RMSD　& RMSF plot--------------------------------------------------#


#plot RMSD & RMSF compare figure
bash /home/duanqi/dq/scripts/rmsd_rmsf_plot.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/rmsd_rmsf_plot.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns

###ACF###------------------------------(RMSD)

# plot ACF of pdb and af RMSD
bash /home/duanqi/dq/scripts/ACF_RMSD.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/ACF_RMSD.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns

# plot RMSD ttest and levene test results of different sample intervals(Using xvg files second column to do t-test and levene-test)
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 1BPI af_SS 5_18 pdb_SS 5_18 0.1 1.0 10 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 2.0 20 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 7.9 10 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 10.0 10 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 2.1 20.0 20 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 20.0 40.0 20 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 30.0 50.0 20 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 1BPI af_SS 5_18 pdb_SS 5_18 0.1 1.5 15 500ns
bash /home/duanqi/dq/scripts/rmsd_p_value.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 0.1 10.0 11 500ns


# plot 10 group heat map (t test and levene test) in different sample interval  (define min_interval & max-interval & interval-groups) and group number (10)
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 10 2.0 10.0 9 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 10.0 10.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 2 1.0 1.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 1.0 1.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 2 7.0 7.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 7.0 7.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 50.0 50.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 45.0 45.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 40.0 40.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 5001 2 55.0 55.0 1 500ns
bash /home/duanqi/dq/scripts/rmsd_sample_heatmap.sh 1BPI af_SS 5_18 pdb_SS 5_18 5001 10 0.7 0.7 1 500ns



#----------------------------------------Energy---------------------------------------------------#


bash /home/duanqi/dq/scripts/energy_analysis.sh 1BPI af_SS 5_18 pdb_SS 5_18
bash /home/duanqi/dq/scripts/energy_analysis.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7

#----------------------------------------FEP---------------------------------------------------#
#plot FEP of PC1 and PC2
bash /home/duanqi/dq/scripts/FEP_PC1_PC2.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns
bash /home/duanqi/dq/scripts/FEP_PC1_PC2.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns


#--------------------------------------animation---------------------------------------------------#


# generate movie of af and pdb trajectory.xtc
bash /home/duanqi/dq/scripts/make_movie.sh 1BPI af_SS 5_18 pdb_SS 5_18 500ns 5 0
bash /home/duanqi/dq/scripts/make_movie.sh 6ZD3 af_triclinic 6_7 pdb_triclinic 6_7 500ns 5 0


#--------------------------------------visualization---------------------------------------------------#
# generate visualization of af and pdb initial structure(need to create analysis directory, having pdb and af structure in results_dir_gromacs)
bash /home/duanqi/dq/scripts/visual_structure.sh 6ZD3 af_triclinic 6_7 0 pdb_triclinic 6_7 5


#-------------------------------------Secondary structure---------------------------------------------------#
# calculate secondary structure of af and pdb initial structure
# bash /home/duanqi/dq/scripts/secondary_structure.sh 1BPI af_SS 5_18 pdb_SS 5_18