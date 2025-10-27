#!/bin/bash

num_procs=$1

for power in 16 18 20 22 24 26 28; do
    for sort_level in 0 1 2 3; do
        num_elements=$((2 ** $power))
        sbatch mpi.grace_job $num_procs $num_elements $sort_level 
    done
done

