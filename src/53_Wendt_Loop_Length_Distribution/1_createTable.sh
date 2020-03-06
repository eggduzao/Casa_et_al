#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CTD
#SBATCH --output=CTD.%A_%a.out
#SBATCH --error=CTD.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-1

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1361560

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/input/1_ctd.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/1_createTable.py $parameters


