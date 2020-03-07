#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RPS
#SBATCH --output=RPS.%A_%a.out
#SBATCH --error=RPS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --array=1-15

# Commands
# sbatch 3_plotSubtraction.sh
# squeue -u egadegu
# scancel 3319963

# Input
inputFileName1="/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/31_rps.txt"
parameters1=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName1`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/3_createFullTable.py $parameters1

# Input
inputFileName2="/usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/input/32_rps.txt"
parameters2=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName2`

# Creating matrix
Rscript /usr/users/egadegu/Projects/Wendt_Stag/Code/45_FigS8B_SubMaps_RaoNeg/3_hic.R $parameters2


