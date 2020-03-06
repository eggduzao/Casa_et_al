#!/bin/bash

# Slurm Parameters
#SBATCH -p fat
#SBATCH --job-name=JUI
#SBATCH --output=JUI.%A_%a.out
#SBATCH --error=JUI.%A_%a.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=10
#SBATCH --mem=256gb
#SBATCH --time=48:00:00
#SBATCH --array=1-6

# Commands
# sbatch 1_applyJuicer.sh
# squeue -u egadegu
# scancel xxxxxxxxxxxxxxx

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/input/1_jui.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/19_Process_Losada_Hic/1_applyJuicer.py $parameters


