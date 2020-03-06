#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CBF
#SBATCH --output=CBF.%A_%a.out
#SBATCH --error=CBF.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-16

# 2MB
# 1-16 1_cbf1.txt # OK 0 err
# 1-1000 1_cbf2.txt # OK 232 err
# 1-1000 1_cbf3.txt # OK 57 err
# 1-992 1_cbf4.txt # OK 0 err

# Commands
# sbatch 1_binMatrixAtRegion.sh
# squeue -u egadegu
# scancel 1914036, 1914121, 1915152, 1916162, 1924606, 1950649

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/1_cbf.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/1_binMatrixAtRegion.py $parameters


