#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CIC
#SBATCH --output=CIC.%A_%a.out
#SBATCH --error=CIC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=5:00:00
#SBATCH --array=1-16

# Commands
# sbatch 1_createTable.sh
# squeue -u egadegu
# scancel 1512157

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/33_ChIP_Intensity_Rank_Comparison/input/1_cic.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/33_ChIP_Intensity_Rank_Comparison/1_createTable.py $parameters


