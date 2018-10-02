#!/bin/bash

for model in {1..7}

  if [ $model -eq 2 ]
  then
    nCPU=3
  else
    nCPU=1
  fi
  
    for dataset in {1..3}
    
      if [ $dataset -eq 1 ]
      then
        MEMORY_LIMIT=20480
        PARTITION="cloud"
        WALLTIME=5
      fi
      
      if [ $dataset -eq 2 ]
      then
        MEMORY_LIMIT=51200
        PARTITION="cloud"
        
        if [ $model -eq 6 ]
        then
          WALLTIME=168
        else
          WALLTIME=24
        fi
        
      fi
      
      if [ $dataset -eq 3 ]
      then
        MEMORY_LIMIT=225280
        PARTITION="physical"
        
        if [ $model -eq 1 ]
        then
          WALLTIME=24
        elif [ $model -eq 2 ]
        then
          WALLTIME=336
        elif [ $model -eq 3 ]
        then
          WALLTIME=336
        elif [ $model -eq 4 ]
        then
          WALLTIME=48
        elif [ $model -eq 5 ]
        then
          WALLTIME=48
        elif [ $model -eq 6 ]
        then
          WALLTIME=336
        elif [ $model -eq 7 ]
        then
          WALLTIME=24
        fi  
      fi
      
        for fold in {1..5}
        
          sbatch /data/cephfs/punim0200/scripts/slurm/master_slurm.slurm $model $dataset $fold $MEMORY_LIMIT $PARTITION $WALLTIME $nCPU
        