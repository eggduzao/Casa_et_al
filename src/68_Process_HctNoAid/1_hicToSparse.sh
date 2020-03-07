#!/bin/bash

# Slurm Parameters
#SBATCH -p fat-fmz
#SBATCH --job-name=HTS
#SBATCH --output=HTS.%A_%a.out
#SBATCH --error=HTS.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200gb
#SBATCH --time=24:00:00
#SBATCH --array=1-14

# Commands
# sbatch 1_hicToSparse.sh
# squeue -u egadegu
# scancel 3229513

# Making temp dir local
export TMPDIR=$TMP_LOCAL

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/input/1_hts.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
time python /usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/1_hicToSparse.py $parameters

