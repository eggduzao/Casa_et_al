#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=EIG
#SBATCH --output=EIG.%A_%a.out
#SBATCH --error=EIG.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-184

# Commands
# sbatch 1_eigenvector.sh
# squeue -u egadegu
# scancel 1428525

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/3_AB_Compartments/input/1_eig.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
juicertools eigenvector -p $parameters


