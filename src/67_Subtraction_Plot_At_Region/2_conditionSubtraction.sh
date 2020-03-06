#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RSC
#SBATCH --output=RSC.%A_%a.out
#SBATCH --error=RSC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-16

# 2MB
# 1-8 2_rcs1.txt # OK
# 1-1000 2_rcs2.txt # OK 259 err
# 1-496 2_rcs3.txt # OK 69 err

# Commands
# sbatch 2_conditionSubtraction.sh
# squeue -u egadegu
# scancel 1914062, 1917184, 1918193, 1924641, 1950849

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/2_rcs.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/2_conditionSubtraction.py $parameters


