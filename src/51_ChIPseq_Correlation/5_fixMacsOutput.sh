#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FMO
#SBATCH --output=FMO.%A_%a.out
#SBATCH --error=FMO.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-2

# Commands
# sbatch 5_fixMacsOutput.sh
# squeue -u egadegu
# scancel 1360996

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/input/5_fmo.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/22_ChIPseq_Correlation/5_fixMacsOutput.py $parameters


