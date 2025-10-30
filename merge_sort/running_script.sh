#!/bin/bash

num_procs=$1

for power in 16 18 20 22 24 26 28; do
    for sort_level in 0 1 2 3; do
        for use_float in 0 1; do
            if [ $use_float -eq 0 ]; then
                dtype="int"
            else
                dtype="float"
            fi
            # Select job file based on num_procs
            case "$num_procs" in
                2)    jobfile="grace-files/mpi-2.grace_job" ;;
                4)    jobfile="grace-files/mpi-4.grace_job" ;;
                8)    jobfile="grace-files/mpi-8.grace_job" ;;
                16)   jobfile="grace-files/mpi-16.grace_job" ;;
                32)   jobfile="grace-files/mpi-32.grace_job" ;;
                64)   jobfile="grace-files/mpi-64.grace_job" ;;
                128)  jobfile="grace-files/mpi-128.grace_job" ;;
                256)  jobfile="grace-files/mpi-256.grace_job" ;;
                512)  jobfile="grace-files/mpi-256.grace_job" ;;
                1024) jobfile="grace-files/mpi-1024.grace_job" ;;
                *)    echo "No job file for $num_procs processes. Skipping."; continue ;;
            esac
            echo "Submitting: sbatch $jobfile $num_procs $power $sort_level $use_float ($dtype)"
            sbatch $jobfile $num_procs $power $sort_level $use_float
            echo "Sleeping for 15 minutes to avoid overloading the scheduler..."
            sleep 900
        done
    done
done