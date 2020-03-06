#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=ALG
#SBATCH --output=ALG.%A_%a.out
#SBATCH --error=ALG.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=48gb
#SBATCH --time=10:00:00
#SBATCH --array=1-20

# Commands
# sbatch 1_alignment.sh
# squeue -u egadegu
# scancel 1250543

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/input/1_alg.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/18_Process_Losada_ChIP/1_alignment.py $parameters


