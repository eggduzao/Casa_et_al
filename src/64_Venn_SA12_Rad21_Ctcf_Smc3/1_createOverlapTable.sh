#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=COT
#SBATCH --output=COT.%A_%a.out
#SBATCH --error=COT.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --time=1:00:00
#SBATCH --array=2-2

# Commands
# sbatch 1_createOverlapTable.sh
# squeue -u egadegu
# scancel 1812152, 1812522

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/35_Venn_SA12_Rad21_Ctcf_Smc3/input/1_cot.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /usr/users/egadegu/Projects/Wendt_Stag/Code/35_Venn_SA12_Rad21_Ctcf_Smc3/1_createOverlapTable.py $parameters


