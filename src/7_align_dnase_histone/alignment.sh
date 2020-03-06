#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=ALIGN
#SBATCH --output=ALIGN.%A_%a.out
#SBATCH --error=ALIGN.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-76

# Commands
# sbatch alignment.sh
# squeue -u egusmao
# scancel 10537566

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
export PATH=$PATH:"/projects/ag-papan/install/BEDTools_2.17.0/bedtools-2.17.0/bin/"
export PATH=$PATH:"/projects/ag-papan/install/juicertools-1.7.6/bin/"
export PYTHONPATH=$PYTHONPATH:"/home/egusmao/.local/lib/python2.7/site-packages"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/"

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Previous_Code/7_align_dnase_histone/rao_alignment/"${SLURM_ARRAY_TASK_ID}".txt"
alignType=`sed '1q;d' $inputFileName`
minQuality=`sed '2q;d' $inputFileName`
ncores=`sed '3q;d' $inputFileName`
fastqGzFileName=`sed '4q;d' $inputFileName`
indexFileName=`sed '5q;d' $inputFileName`
tempLocation=`sed '6q;d' $inputFileName`
outputLocation=`sed '7q;d' $inputFileName`

# Creating matrix
time python /projects/ag-papan/eduardo/Papantonis_Stag/Previous_Code/7_align_dnase_histone/alignment.py $alignType $minQuality $ncores $fastqGzFileName $indexFileName $tempLocation $outputLocation 


