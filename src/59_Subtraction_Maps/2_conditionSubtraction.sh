#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RSM
#SBATCH --output=RSM.%A_%a.out
#SBATCH --error=RSM.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=10:00:00
#SBATCH --array=1-176

# Commands
# sbatch 2_conditionSubtraction.sh
# squeue -u egadegu
# scancel 1810232

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/2_rsm.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/2_conditionSubtraction.py $parameters


