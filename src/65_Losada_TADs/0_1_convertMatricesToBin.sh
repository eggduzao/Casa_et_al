#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CMB
#SBATCH --output=CMB.%A_%a.out
#SBATCH --error=CMB.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=1-69

# Commands
# sbatch 0_1_convertMatricesToBin.sh
# squeue -u egadegu
# scancel 1696253

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/input/01_cmb.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Execution
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/0_1_convertMatricesToBin.py $parameters


