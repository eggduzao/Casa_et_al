#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=merge_hic_3B9
#SBATCH --output=merge_hic_3B9.out
#SBATCH --error=merge_hic_3B9.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=42gb
#SBATCH --time=50:00:00
#SBATCH --account=UniKoeln

# Commands
# sbatch merge.sh
# squeue -u njosipov
# scancel xxxxxxx

# Packages
module load jdk/1.8.0_60

# Input
inFileName1="/projects/ag-papan/Yulia/HiC-STAG1-2/79644/aligned/new_merged_nodups.txt"
inFileName2="/projects/ag-papan/Yulia/HiC-STAG1-2-R1/75728/aligned/merged_nodups.txt"

# Test
cd /projects/ag-papan/eduardo/Papantonis_Stag/3B9_5+/
/projects/ag-papan/eduardo/Papantonis_Stag/3B9_5+/mega.sh $inFileName1 $inFileName2


