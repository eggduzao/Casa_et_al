#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CBF
#SBATCH --output=CBF.%A_%a.out
#SBATCH --error=CBF.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --array=1-9

# Commands
# sbatch 1_binMatrixAtRegion.sh
# squeue -u egadegu
# scancel 3319133

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/1_cbf.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/1_binMatrixAtRegion.py $parameters


