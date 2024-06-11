#!/bin/bash

# Initialize modules
source /etc/profile

# Loading modules
module purge
module load julia/1.9.2
module load gurobi/gurobi-1000

echo "My SLURM_ARRAY_TASK_ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

julia script.jl $LLSUB_RANK $LLSUB_SIZE