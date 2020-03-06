#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RPE
#SBATCH --output=RPE.%A_%a.out
#SBATCH --error=RPE.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-8

# Commands
# sbatch 2_readsPerExon.sh
# squeue -u egadegu
# scancel 1741212

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/input/2_rpe.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/37_Basic_Statistics/2_readsPerExon.py $parameters


