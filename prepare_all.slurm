#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=debug
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 16000
#SBATCH --time 20:00:00
#SBATCH --array=0-999
#SBATCH --output=exelogs/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/_masif_precompute.%A_%a.err

i=0

while read p; do

  if [ $(( i % 1000 )) == ${SLURM_ARRAY_TASK_ID} ]; then
    FIELD1=$(echo $p)
    start_time=$(date +%s)
    python prepare_all.py $FIELD1
    end_time=$(date +%s)
    cost_time=$[ $end_time-$start_time ]
    #echo "$i	${SLURM_ARRAY_TASK_ID}  $p	$cost_time"
  fi
  i=$((i+1))
#done < data/pdbid_chain/RBP03.txt
done < data/pdbid_chain/RBP407_09.txt
