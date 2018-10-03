#!/bin/bash
#
# Set compute node type. Read in PARTITION command line argument ($5)
#SBATCH -p $5
#
# Set required nodes.
#SBATCH --nodes=1
#
# Set up nodes for job. Read in nCPU command line argument ($7)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$7
#
# Set memory limits.  Read in MEMORY_LIMIT command line argument ($4)
#SBATCH --mem=$4
#
# Set walltime. Read in WALLTIME command line argument ($6). Set in Hours:Minutes:Seconds format
#SBATCH --time=$6:00:00
#
# Set up job name. Format of JPred_ModelID_DatasetID_FoldID
#SBATCH --job-name=JPred_$1_$2_$3
#
# Set up output file name. Format of JPred_ModelID_DatasetID_FoldID.out
#SBATCH --output=outputs/slurm_outputs/JPred_$1_$2_$3.out
#
# Set email notifications
#SBATCH --mail-user=davidpw@student.unimelb.edu.au
#SBATCH --mail-type=ALL,TIME_LIMIT_50,TIME_LIMIT_80,TIME_LIMIT_90

model=$1
dataset=$2
fold=$3

module purge

module load PROJ/5.0.0-intel-2017.u2
module load GDAL/2.2.3-intel-2017.u2
module load JAGS/4.3.0-intel-2017.u2
module load R/3.5.0-intel-2017.u2-X11-20171023

cd /data/cephfs/punim0200/JSDM_Prediction

Rscript --vanilla scripts/master/JSDM_Prediction_Master.R $model $dataset $fold