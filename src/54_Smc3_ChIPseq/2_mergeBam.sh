#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MEB
#SBATCH --output=MEB.%A_%a.out
#SBATCH --error=MEB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --array=1-2

# Commands
# sbatch 2_mergeBam.sh
# squeue -u egadegu
# scancel 1367557

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/2_meb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/2_mergeBam.py $parameters


