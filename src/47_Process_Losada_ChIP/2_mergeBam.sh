#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MEB
#SBATCH --output=MEB.%A_%a.out
#SBATCH --error=MEB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=5:00:00
#SBATCH --array=1-14

# Commands
# sbatch 2_mergeBam.sh
# squeue -u egadegu
# scancel 1257508

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/input/2_meb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/2_mergeBam.py $parameters


