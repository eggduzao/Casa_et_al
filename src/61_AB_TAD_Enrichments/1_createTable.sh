#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=ATEB
#SBATCH --output=ATEB.%A_%a.out
#SBATCH --error=ATEB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=5:00:00
#SBATCH --array=1-12

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1441967, 1447259

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/32_AB_TAD_Enrichments/input/1_ateb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/32_AB_TAD_Enrichments/1_createTable.py $parameters


