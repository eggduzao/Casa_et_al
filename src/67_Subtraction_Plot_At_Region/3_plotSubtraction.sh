#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RPS
#SBATCH --output=RPS.%A_%a.out
#SBATCH --error=RPS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=2:00:00
#SBATCH --array=1-176

# 2MB
# 1-8 31_rps1  # OK
# 1-1000 31_rps2 # OK
# 1-496 31_rps3 # OK

# 2MB
# 1-8 32_rps1  # OK
# 1-1000 32_rps2 # OK
# 1-496 32_rps3 # OK

# Commands
# sbatch 3_plotSubtraction.sh
# squeue -u egadegu
# scancel 1914096, 1914112, 1918700, 1919710, 1920207, 1921211, 1924684, 1950905, 1951481

# Input
#inputFileName1="/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/31_rps.txt"
#parameters1=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName1`

# Creating matrix
#python /usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/3_createFullTable.py $parameters1

# Input
inputFileName2="/usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/input/32_rps.txt"
parameters2=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName2`

# Creating matrix
Rscript /usr/users/egadegu/Projects/Wendt_Stag/Code/38_Subtraction_Plot_At_Region/3_hic.R $parameters2


