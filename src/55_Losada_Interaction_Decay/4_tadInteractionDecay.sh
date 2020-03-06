#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=TID
#SBATCH --output=TID.%A_%a.out
#SBATCH --error=TID.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-368

# Commands
# sbatch 4_tadInteractionDecay.sh
# squeue -u egadegu
# scancel 1385784

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/4_tid.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/4_tadInteractionDecay.py $parameters


