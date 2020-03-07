#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RCS
#SBATCH --output=RCS.%A_%a.out
#SBATCH --error=RCS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-6

# Commands
# sbatch 2_conditionSubtraction.sh
# squeue -u egadegu
# scancel 3319504

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/2_rcs.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/2_conditionSubtraction.py $parameters


