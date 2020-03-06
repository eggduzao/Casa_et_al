#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CBW
#SBATCH --output=CBW.%A_%a.out
#SBATCH --error=CBW.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --array=1-2

# Commands
# sbatch 6_createBigWig.sh
# squeue -u egadegu
# scancel 1368046

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/6_cbw.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/6_createBigWig.py $parameters


