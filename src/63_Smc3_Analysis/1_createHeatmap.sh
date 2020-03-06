#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CHS
#SBATCH --output=CHS.%A_%a.out
#SBATCH --error=CHS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-4

# Commands
# sbatch 1_createHeatmap.sh
# squeue -u egadegu
# scancel 1551119, 1818972, 1943059

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/34_Smc3_Analysis/input/1_chs.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/34_Smc3_Analysis/1_createHeatmap.py $parameters


