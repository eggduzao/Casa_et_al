#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FSP2
#SBATCH --output=FSP2.%A_%a.out
#SBATCH --error=FSP2.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=3:00:00
#SBATCH --array=1001-2000

# 1 31_fsp
# 1-1000 31_fsp1  # OK
# 1001-2000 31_fsp1 # OK
# 1-1000 31_fsp2 # OK
# 1001-1036 31_fsp2 # OK

# 2 32_fsp
# 1-1000 32_fsp1  # OK
# 1001-2000 32_fsp1 # RUNNING
# 1-1000 32_fsp2 # 
# 1001-1036 32_fsp2 # 

# Commands
# sbatch 3_createFullTable.sh
# squeue -u egadegu
# scancel 1811085, 1812128, 1814822, 1818492, 1821881, 1850745, 1852815, 1888286, 

# Input
#inputFileName1="/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/31_fsp2.txt"
#parameters1=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName1`

# Creating matrix
#python /usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/3_createFullTable.py $parameters1

# Input
inputFileName2="/usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/input/32_fsp1.txt"
parameters2=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName2`

# Creating matrix
Rscript /usr/users/egadegu/Projects/Wendt_Stag/Code/30_Subtraction_Maps/3_hic.R $parameters2


