#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RPP
#SBATCH --output=RPP.%A_%a.out
#SBATCH --error=RPP.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-11

# Commands
# sbatch 1_readsPerPeak.sh
# squeue -u egadegu
# scancel 1741211

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/input/1_rpp.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/1_readsPerPeak.py $parameters


