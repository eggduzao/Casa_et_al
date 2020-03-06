#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CTT
#SBATCH --output=CTT.%A_%a.out
#SBATCH --error=CTT.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-1

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1390977

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/27_Number_TADs_per_Chromosome/input/1_ctt.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/27_Number_TADs_per_Chromosome/1_createTable.py $parameters


