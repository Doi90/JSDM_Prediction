#!/bin/bash

for model in {1..7}
do
  if [ $model = "2" ]
  then
    nCPU=3
  else
    nCPU=1
  fi
  
    for dataset in {1..20}
    do
	    MEMORY_LIMIT=20480
      PARTITION="cloud,physical"
      WALLTIME=7	
      
        for fold in {1..5}
        do
		      WALLTIME_FORMAT="$WALLTIME-00"
		      JOB_NAME="S_${model}_${dataset}_${fold}"
		      OUT_NAME="/data/cephfs/punim0200/JSDM_Prediction/outputs/slurm_outputs/S_${model}_${dataset}_${fold}.out"
		  
          sbatch -p $PARTITION --cpus-per-task=$nCPU --mem=$MEMORY_LIMIT --time=$WALLTIME_FORMAT --job-name=$JOB_NAME --output=$OUT_NAME /data/cephfs/punim0200/JSDM_Prediction/scripts/slurm/master_slurm_simulation.slurm $model $dataset $fold 
        done
	done	  
done        