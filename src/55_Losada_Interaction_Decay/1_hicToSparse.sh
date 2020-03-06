#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=HTS
#SBATCH --output=HTS.%A_%a.out
#SBATCH --error=HTS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=5:00:00
#SBATCH --array=1-24

# Commands
# sbatch 1_hicToSparse.sh
# squeue -u egadegu
# scancel 1381845

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/input/1_hts.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/26_Losada_Interaction_Decay/1_hicToSparse.py $parameters


