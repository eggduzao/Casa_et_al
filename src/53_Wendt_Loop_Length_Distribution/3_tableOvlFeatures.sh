#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CTC
#SBATCH --output=CTC.%A_%a.out
#SBATCH --error=CTC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=10:00:00
#SBATCH --array=1-1

# Commands
# sbatch 3_tableOvlFeatures.sh
# squeue -u egadegu
# scancel 1714579

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/input/3_ctc.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/24_Wendt_Loop_Length_Distribution/3_tableOvlFeatures.py $parameters


