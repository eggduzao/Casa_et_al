#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=EIG
#SBATCH --output=EIG.%A_%a.out
#SBATCH --error=EIG.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48gb
#SBATCH --time=3:00:00
#SBATCH --array=1-46

# Commands
# sbatch 1_eigenvector.sh
# squeue -u egadegu
# scancel 3262805, 3264627

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/1_eig.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
juicertools eigenvector -p $parameters


