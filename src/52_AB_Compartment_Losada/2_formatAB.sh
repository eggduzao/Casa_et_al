#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=FAB
#SBATCH --output=FAB.%A_%a.out
#SBATCH --error=FAB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --array=1-69

# Commands
# sbatch 2_formatAB.sh
# squeue -u egadegu
# scancel 1361692

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/23_AB_Compartment_Losada/input/2_fab.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/23_AB_Compartment_Losada/2_formatAB.py $parameters


