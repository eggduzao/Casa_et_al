#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CPT
#SBATCH --output=CPT.%A_%a.out
#SBATCH --error=CPT.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48gb
#SBATCH --time=2:00:00
#SBATCH --array=1-8

# Commands
# sbatch 3_createPlotTable.sh
# squeue -u egadegu
# scancel 3265089

# Input 1
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/31_cpt.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating table
python /usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/3_createPlotTable.py $parameters

# Input 2
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/input/32_cpt.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Plotting table
Rscript /usr/users/egadegu/Projects/Wendt_Stag/Code/41_FigS9A_AB_RaoNeg/3_createPlotTable.R $parameters

