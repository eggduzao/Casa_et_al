#!/bin/bash

# Slurm Parameters
#SBATCH -p fat
#SBATCH --job-name=RMH
#SBATCH --output=RMH.%A_%a.out
#SBATCH --error=RMH.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256gb
#SBATCH --time=24:00:00
#SBATCH --array=1-3

# Commands
# sbatch 2_mergeHic.sh
# squeue -u egadegu
# scancel 1257545

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/input/2_rmh.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/2_mergeHic.py $parameters


