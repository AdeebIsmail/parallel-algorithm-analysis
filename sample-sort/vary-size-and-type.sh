#!/bin/bash

num_procs=$1

for data_type in int double; do
  for power in 16 18 20 22 24 26 28; do
    for sort_level in random sorted reversed perturbed1; do
      num_elements=$((2 ** $power))
      sbatch ./grace-batch-files/${num_procs}-mpi.grace_job $num_elements $data_type $sort_level "weak"
    done
  done
done
