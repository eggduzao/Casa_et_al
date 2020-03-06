#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MLG
#SBATCH --output=MLG.%A_%a.out
#SBATCH --error=MLG.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-4

# Commands
# sbatch 1_metaLoopOGrams.sh
# squeue -u egusmao
# scancel 10674328

# Modules
module add python/2.7.5-V2
module add intel/17.0_gnu_5.1.0
module add bowtie2/2.2.9
module add samtools/1.6
module add R/3.3.3_intel_mkl
export PATH=$PATH:"/home/egusmao/.local/bin"
export PATH=$PATH:"/projects/ag-papan/install/FastQC/"
export PATH=$PATH:"/projects/ag-papan/install/cutadapt-1.15/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/TrimGalore-0.4.3/"
export PATH=$PATH:"/projects/ag-papan/install/sratoolkit.2.8.2-1-ubuntu64/bin/"
export PATH=$PATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/bin/"
export PATH=$PATH:"/projects/ag-papan/install/bedtools-2.25.0/bin/"
export PATH=$PATH:"/projects/ag-papan/install/juicertools-1.7.6/bin/"
export PATH=$PATH:"/projects/ag-papan/install/chromhmm-1.15/bin/"
export PATH=$PATH:"/home/egusmao/software/ibm-cbc-genomic-tools-master/bin/"
export PYTHONPATH=$PYTHONPATH:"/home/egusmao/.local/lib/python2.7/site-packages"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/"

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/8_Meta_Plots/input/"${SLURM_ARRAY_TASK_ID}"_mlg.txt"
loopBins=`sed '1q;d' $inputFileName`
resolution=`sed '2q;d' $inputFileName`
regionsFile1Name=`sed '3q;d' $inputFileName`
regionsFile2Name=`sed '4q;d' $inputFileName`
hicFile1Name=`sed '5q;d' $inputFileName`
hicFile2Name=`sed '6q;d' $inputFileName`
outputFileName=`sed '7q;d' $inputFileName`

# Execution
python /projects/ag-papan/eduardo/Papantonis_Stag/Code/8_Meta_Plots/1_metaLoopOGrams.py $loopBins $resolution $regionsFile1Name $regionsFile2Name $hicFile1Name $hicFile2Name $outputFileName


