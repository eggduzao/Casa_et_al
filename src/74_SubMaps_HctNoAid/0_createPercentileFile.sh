#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CPF
#SBATCH --output=CPF.%A_%a.out
#SBATCH --error=CPF.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=60gb
#SBATCH --time=5:00:00
#SBATCH --array=1-1

# Commands
# sbatch 0_createPercentileFile.sh
# squeue -u egadegu
# scancel 3313554

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/0_cpf.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/0_createPercentileFile.py $parameters


