#!/bin/bash

# Slurm Parameters
#SBATCH -p fat-fmz
#SBATCH --job-name=RHI
#SBATCH --output=RHI.%A_%a.out
#SBATCH --error=RHI.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200gb
#SBATCH --time=12:00:00
#SBATCH --array=1-322

# Commands
# sbatch 3_chromToHic.sh
# squeue -u egadegu
# scancel 3254208

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/input/3_rhi.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/3_chromToHic.py $parameters


