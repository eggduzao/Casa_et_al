#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FAB
#SBATCH --output=FAB.%A_%a.out
#SBATCH --error=FAB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48gb
#SBATCH --time=2:00:00
#SBATCH --array=1-46

# Commands
# sbatch 2_formatAB.sh
# squeue -u egadegu
# scancel 3263998, 3264948

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/2_fab.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/2_formatAB.py $parameters


