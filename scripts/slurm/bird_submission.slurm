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
start_sample=$4
end_sample=$5

module purge

module load PROJ/5.0.0-GCC-6.2.0
module load GDAL/2.2.3-GCC-6.2.0
module load JAGS/4.3.0-GCC-6.2.0
module load R/3.5.0-GCC-6.2.0

cd /data/cephfs/punim0200/JSDM_Prediction

Rscript --vanilla scripts/master/JSDM_Bird_Prediction_Master.R $model $dataset $fold $start_sample $end_sample