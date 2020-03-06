#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=MAC
#SBATCH --output=MAC.%A_%a.out
#SBATCH --error=MAC.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=11

# Commands
# sbatch 4_macs2.sh
# squeue -u egusmao
# scancel 10873186, 10931925

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
inputFileName="/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/4_mac.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/4_macs2.py $parameters


