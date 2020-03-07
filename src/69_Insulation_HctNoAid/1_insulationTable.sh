#!/bin/bash

# Slurm Parameters
#SBATCH -p fat-fas
#SBATCH --job-name=CIP
#SBATCH --output=CIP.%A_%a.out
#SBATCH --error=CIP.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=200gb
#SBATCH --time=24:00:00
#SBATCH --array=1-14

# Commands
# sbatch 1_insulationTable.sh
# squeue -u egadegu
# scancel 3266152

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/40_Fig4G_Insulation_RaoNeg/input/1_cip.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/40_Fig4G_Insulation_RaoNeg/1_insulationTable.py $parameters


