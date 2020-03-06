#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MAC
#SBATCH --output=MAC.%A_%a.out
#SBATCH --error=MAC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=3:00:00
#SBATCH --array=1-2

# Commands
# sbatch 4_macs2.sh
# squeue -u egadegu
# scancel 1360900

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/4_mac.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/4_macs2.py $parameters


