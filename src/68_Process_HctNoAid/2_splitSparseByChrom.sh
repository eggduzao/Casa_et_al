#!/bin/bash

# Slurm Parameters
#SBATCH -p medium-fmz
#SBATCH --job-name=SPL
#SBATCH --output=SPL.%A_%a.out
#SBATCH --error=SPL.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --array=1-322

# Commands
# sbatch 2_splitSparseByChrom.sh
# squeue -u egadegu
# scancel 3244255

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/input/2_spl.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/2_splitSparseByChrom.py $parameters


