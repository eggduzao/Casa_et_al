#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=merge_hic_69
#SBATCH --output=merge_hic_69.out
#SBATCH --error=merge_hic_69.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=42gb
#SBATCH --time=50:00:00
#SBATCH --account=UniKoeln

# Commands
# sbatch merge.sh
# squeue -u njosipov
# scancel 9332036

# Packages
module load jdk/1.8.0_60

# Input
inFileName1="/projects/ag-papan/Yulia/HiC-STAG1-2/79646/aligned/merged_nodups.txt"
inFileName2="/projects/ag-papan/Yulia/HiC-STAG1-2-R1/75726/aligned/merged_nodups.txt"

# Test
cd /projects/ag-papan/eduardo/Papantonis_Stag/69_127+/
/projects/ag-papan/eduardo/Papantonis_Stag/69_127+/mega.sh $inFileName1 $inFileName2


