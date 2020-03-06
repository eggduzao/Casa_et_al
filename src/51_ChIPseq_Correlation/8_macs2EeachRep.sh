#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MACR
#SBATCH --output=MACR.%A_%a.out
#SBATCH --error=MACR.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=3:00:00
#SBATCH --array=3-4

# Commands
# sbatch 8_macs2EeachRep.sh
# squeue -u egadegu
# scancel 1811789

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/8_macr.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/8_macs2EeachRep.py $parameters


