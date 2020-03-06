#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=CBW
#SBATCH --output=CBW.%A_%a.out
#SBATCH --error=CBW.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --time=24:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=15-16

# Commands
# sbatch 6_createBigWig.sh
# squeue -u egusmao
# scancel 10873199, 10931929

# Modules
module add python/2.7.5-V2
module add bowtie2/2.2.9
module add samtools/1.6
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
inputFileName="/projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/input/6_cbw.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Creating matrix
python /projects/ag-papan/eduardo/Wendt_Stag/Code/9_Process_Losada_Data/6_createBigWig.py $parameters


