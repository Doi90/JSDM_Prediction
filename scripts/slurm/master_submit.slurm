#!/bin/bash

for model in {1..7}
do
  if [ $model = "2" ]
  then
    nCPU=3
  else
    nCPU=1
  fi
  
    for dataset in {1..3}
    do
	  if [ $dataset = "1" ]
      then
        MEMORY_LIMIT=20480
        PARTITION="cloud"
        WALLTIME=7	
      fi
      
      if [ $dataset = "2" ]
      then
        MEMORY_LIMIT=20480
        PARTITION="cloud"
        WALLTIME=10
      fi
      
      if [ $dataset = "3" ]
      then
        
        if [ $model = "1" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="cloud"
          WALLTIME=30
        elif [ $model = "2" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="physical"
          WALLTIME=30
        elif [ $model = "3" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="physical"
          WALLTIME=30
        elif [ $model = "4" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="physical"
          WALLTIME=30
        elif [ $model = "5" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="physical"
          WALLTIME=30
        elif [ $model = "6" ]
        then
          MEMORY_LIMIT=102400
          PARTITION="physical"
          WALLTIME=30
        elif [ $model = "7" ]
        then
          MEMORY_LIMIT=51200
          PARTITION="cloud"
          WALLTIME=7
        fi  
      fi
      
        for fold in {1..5}
        do
		  WALLTIME_FORMAT="$WALLTIME-00"
		  JOB_NAME="JP_${model}_${dataset}_${fold}"
		  OUT_NAME="/data/cephfs/punim0200/JSDM_Prediction/outputs/slurm_outputs/JP_${model}_${dataset}_${fold}.out"
		  
          sbatch -p $PARTITION --cpus-per-task=$nCPU --mem=$MEMORY_LIMIT --time=$WALLTIME_FORMAT --job-name=$JOB_NAME --output=$OUT_NAME /data/cephfs/punim0200/JSDM_Prediction/scripts/slurm/master_slurm.slurm $model $dataset $fold 
        done
	done	  
done        