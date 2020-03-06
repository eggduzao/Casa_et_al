#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=WLH
#SBATCH --output=WLH.%A_%a.out
#SBATCH --error=WLH.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --array=1-4

# Commands
# sbatch 1_heatmapTable.sh
# squeue -u egadegu
# scancel 1359353, 1537352

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/21_Comp_Losada_ChIP_Heatmap/input/1_wlh.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/21_Comp_Losada_ChIP_Heatmap/1_heatmapTable.py $parameters


