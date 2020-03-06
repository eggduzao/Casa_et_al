#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=ABTB
#SBATCH --output=ABTB.%A_%a.out
#SBATCH --error=ABTB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=5:00:00
#SBATCH --array=1-2

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1442205, 1447273

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/31_Average_Diff_Interactions_AB/input/1_abtb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/31_Average_Diff_Interactions_AB/1_createTable.py $parameters


