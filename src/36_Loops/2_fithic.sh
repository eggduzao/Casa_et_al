#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=RFH
#SBATCH --output=RFH.%A_%a.out
#SBATCH --error=RFH.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=1:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-368

# Commands
# sbatch 2_fithic.sh
# squeue -u egusmao
# scancel 10581350

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
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/5_Loops/input/"${SLURM_ARRAY_TASK_ID}"_rfh.txt"
resolution=`sed '1q;d' $inputFileName`
lowerBound=`sed '2q;d' $inputFileName`
upperBound=`sed '3q;d' $inputFileName`
equalOcBinNb=`sed '4q;d' $inputFileName`
splineTimes=`sed '5q;d' $inputFileName`
contactFileName=`sed '6q;d' $inputFileName`
fragmentFileName=`sed '7q;d' $inputFileName`
outputName=`sed '8q;d' $inputFileName`
outputLocation=`sed '9q;d' $inputFileName`

# Execution
python /projects/ag-papan/eduardo/Papantonis_Stag/Code/5_Loops/2_fithic.py $resolution $lowerBound $upperBound $equalOcBinNb $splineTimes $contactFileName $fragmentFileName $outputName $outputLocation


