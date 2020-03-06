#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=TID
#SBATCH --output=TID.out
#SBATCH --error=TID.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=5:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-138

# Commands
# sbatch 1_tad_interaction_decay.sh
# squeue -u egusmao
# scancel 9942350

# Modules
module add python/2.7.5-V2
module add intel/17.0_gnu_5.1.0
module add bowtie2/2.2.9
module add samtools/1.6
export PATH=$PATH:"/home/egusmao/.local/bin"
export PATH=$PATH:"/projects/ag-papan/install/FastQC/"
export PATH=$PATH:"/projects/ag-papan/install/cutadapt-1.15/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/TrimGalore-0.4.3/"
export PATH=$PATH:"/projects/ag-papan/install/sratoolkit.2.8.2-1-ubuntu64/bin/"
export PATH=$PATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/bedtools-2.25.0/bin/"
export PATH=$PATH:"/projects/ag-papan/install/juicertools-1.7.6/bin/"
export PYTHONPATH=$PYTHONPATH:"/home/egusmao/.local/lib/python2.7/site-packages"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/"

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/22_tad_interaction_decay/1_input/"${SLURM_ARRAY_TASK_ID}".txt"
chromosome=`sed '1q;d' $inputFileName`
resolution=`sed '2q;d' $inputFileName`
chromSizesFileName=`sed '3q;d' $inputFileName`
matrixMinusFileName=`sed '4q;d' $inputFileName`
matrixPlusFileName=`sed '5q;d' $inputFileName`
outputFileName=`sed '6q;d' $inputFileName`

# Execution
python /projects/ag-papan/eduardo/Papantonis_Stag/Code/22_tad_interaction_decay/1_tad_interaction_decay.py $chromosome $resolution $chromSizesFileName $matrixMinusFileName $matrixPlusFileName $outputFileName


