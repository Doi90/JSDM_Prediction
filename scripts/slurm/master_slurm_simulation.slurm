#!/bin/bash
#
# Set required nodes.
#SBATCH --nodes=1
#
# Set up nodes for job. 
#SBATCH --ntasks=1
#
# Set email notifications
#SBATCH --mail-user=davidpw@student.unimelb.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90

model=$1
dataset=$2
fold=$3

module purge

module load PROJ/5.0.0-GCC-6.2.0
module load GDAL/2.2.3-GCC-6.2.0
module load JAGS/4.3.0-spartan_gcc-6.2.0
module load R/3.5.2-spartan_gcc-6.2.0-bare

cd /data/cephfs/punim0200/JSDM_Prediction

Rscript --vanilla scripts/master/JSDM_Prediction_Simulation_Master.R $model $dataset $fold