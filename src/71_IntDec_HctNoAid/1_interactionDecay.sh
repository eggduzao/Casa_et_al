#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CID
#SBATCH --output=CID.%A_%a.out
#SBATCH --error=CID.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-46

# Commands
# sbatch 1_interactionDecay.sh
# squeue -u egadegu
# scancel 3280596, 3281384

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/43_FigS9B_IntDec_RaoNeg/input/1_cid.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/43_FigS9B_IntDec_RaoNeg/1_interactionDecay.py $parameters


