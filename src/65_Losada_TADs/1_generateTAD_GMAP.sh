#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=GMAP
#SBATCH --output=GMAP.%A_%a.out
#SBATCH --error=GMAP.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=5:00:00
#SBATCH --array=1-69

# Commands
# sbatch 1_generateTAD_GMAP.sh
# squeue -u egadegu
# scancel 1714710

# Input
inputFileName="/usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/input/1_gmap.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Execution
tempDir="/scratch/egadegu/GMAP/"${SLURM_ARRAY_TASK_ID}"/"
mkdir -p $tempDir
cp /usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/1_generateTAD_GMAP.R $tempDir
cd $tempDir
Rscript 1_generateTAD_GMAP.R $parameters


