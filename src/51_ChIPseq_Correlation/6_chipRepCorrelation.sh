#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CRC
#SBATCH --output=CRC.%A_%a.out
#SBATCH --error=CRC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb
#SBATCH --time=10:00:00
#SBATCH --array=1-1

# Commands
# sbatch 6_chipRepCorrelation.sh
# squeue -u egadegu
# scancel 1367416

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/6_crc.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/6_chipRepCorrelation.py $parameters


