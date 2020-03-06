#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CBS
#SBATCH --output=CBS.%A_%a.out
#SBATCH --error=CBS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=5:00:00
#SBATCH --array=1-160

# Commands
# sbatch 1_createBinStdMatrices.sh
# squeue -u egadegu
# scancel 1809846

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/1_cbs.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/1_createBinStdMatrices.py $parameters


