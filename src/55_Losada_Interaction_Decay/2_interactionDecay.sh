#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RID
#SBATCH --output=RID.%A_%a.out
#SBATCH --error=RID.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-552

# Commands
# sbatch 2_interactionDecay.sh
# squeue -u egadegu
# scancel 1384132

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/2_rid.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/2_interactionDecay.py $parameters


