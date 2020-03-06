#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FBW
#SBATCH --output=FBW.%A_%a.out
#SBATCH --error=FBW.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=5:00:00
#SBATCH --array=1-23

# Commands
# sbatch 7_fixBigWig.sh
# squeue -u egadegu
# scancel 1941172

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/7_fbw.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/7_fixBigWig.py $parameters


