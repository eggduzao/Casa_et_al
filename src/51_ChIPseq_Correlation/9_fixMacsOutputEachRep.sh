#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FMOR
#SBATCH --output=FMOR.%A_%a.out
#SBATCH --error=FMOR.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-4

# Commands
# sbatch 9_fixMacsOutputEachRep.sh
# squeue -u egadegu
# scancel 1812137

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/9_fmor.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/9_fixMacsOutputEachRep.py $parameters


