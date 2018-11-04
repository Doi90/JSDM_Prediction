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
        MEMORY_LIMIT=51200
        PARTITION="cloud"
        WALLTIME=3	
      fi
      
      if [ $dataset = "2" ]
      then
        MEMORY_LIMIT=51200
        PARTITION="cloud"
        
        if [ $model = "6" ]
        then
          WALLTIME=10
        else
          WALLTIME=6
        fi
        
      fi
      
      if [ $dataset = "3" ]
      then
        MEMORY_LIMIT=225280
        PARTITION="physical"
        
        if [ $model = "1" ]
        then
          WALLTIME=7
        elif [ $model = "2" ]
        then
          WALLTIME=20
        elif [ $model = "3" ]
        then
          WALLTIME=20
        elif [ $model = "4" ]
        then
          WALLTIME=10
        elif [ $model = "5" ]
        then
          WALLTIME=10
        elif [ $model = "6" ]
        then
          WALLTIME=20
        elif [ $model = "7" ]
        then
          WALLTIME=6
        fi  
      fi
      
        for fold in {1..5}
        do
		  WALLTIME_FORMAT="$WALLTIME-00"
		  JOB_NAME="JPred_$model_$dataset_$fold"
		  OUT_NAME="/data/cephfs/punim0200/JSDM_Prediction/outputs/slurm_outputs/JPred_$model_$dataset_$fold.out"
		  
          sbatch -p $PARTITION --cpus-per-task=$nCPU --mem=$MEMORY_LIMIT --time=$WALLTIME_FORMAT --job-name=$JOB_NAME --output=$OUT_NAME /data/cephfs/punim0200/JSDM_Prediction/scripts/slurm/master_slurm.slurm $model $dataset $fold 
        done
	done	  
done        