#!/bin/bash
#
# Set resource requirements
#SBATCH -p=cloud,physical
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000
#SBATCH --time=168:00:00
#
# Job parameters
#SBATCH --job-name=plot_df
#SBATCH --output=/data/cephfs/punim0200/JSDM_Prediction/slurm/plot_df.out
#
# Set email notifications
#SBATCH --mail-user=davidpw@student.unimelb.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90

module purge

module load PROJ/5.0.0-GCC-6.2.0
module load GDAL/2.2.3-GCC-6.2.0
module load JAGS/4.3.0-spartan_gcc-6.2.0
module load R/3.5.2-spartan_gcc-6.2.0-bare

cd /data/cephfs/punim0200/JSDM_Prediction

Rscript --vanilla scripts/analysis/JSDM_Prediction_circos_figures_chapter3.R
