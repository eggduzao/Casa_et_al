#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CTB
#SBATCH --output=CTB.%A_%a.out
#SBATCH --error=CTB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-4

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1424856

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/28_TAD_Border_Strength/input/1_ctb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/28_TAD_Border_Strength/1_createTable.py $parameters


