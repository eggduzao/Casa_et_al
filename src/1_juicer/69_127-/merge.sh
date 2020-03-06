#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=merge_hic_1
#SBATCH --output=merge_hic_1.out
#SBATCH --error=merge_hic_1.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=42gb
#SBATCH --time=50:00:00
#SBATCH --account=UniKoeln

# Packages
module load jdk/1.8.0_60

# Input
inFileName1="/projects/ag-papan/Yulia/HiC-STAG1-2/79645/aligned/merged_nodups_old.txt"
inFileName2="/projects/ag-papan/Yulia/HiC-STAG1-2-R1/75727/aligned/merged_nodups.txt"

# Test
cd /projects/ag-papan/eduardo/Papantonis_Collaboration/Kargapolova_HicCorrelation/merged_hic_files/69_127-/
/projects/ag-papan/eduardo/Papantonis_Collaboration/Kargapolova_HicCorrelation/merged_hic_files/69_127-/mega.sh $inFileName1 $inFileName2


