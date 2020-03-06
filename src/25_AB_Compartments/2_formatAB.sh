#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=AB
#SBATCH --output=AB.out
#SBATCH --error=AB.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=5:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-184

# Commands
# sbatch 2_formatAB.sh
# squeue -u egusmao
# scancel 10038353

# Modules
module add python/2.7.5-V2
module add intel/17.0_gnu_5.1.0
module add bowtie2/2.2.9
module add samtools/1.6
module add R/3.3.3_intel_mkl
export PATH="/home/egusmao/.local/bin":$PATH
export PATH="/projects/ag-papan/install/FastQC/":$PATH
export PATH="/projects/ag-papan/install/cutadapt-1.15/inst/bin/":$PATH
export PATH="/projects/ag-papan/install/TrimGalore-0.4.3/":$PATH
export PATH="/projects/ag-papan/install/sratoolkit.2.8.2-1-ubuntu64/bin/":$PATH
export PATH="/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/bin/":$PATH
export PATH="/projects/ag-papan/install/bedtools-2.25.0/bin/":$PATH
export PATH="/projects/ag-papan/install/juicertools-1.7.6/bin/":$PATH
export PATH="/projects/ag-papan/install/chromhmm-1.15/bin/":$PATH
export PYTHONPATH="/home/egusmao/.local/lib/python2.7/site-packages":$PYTHONPATH
export PYTHONPATH="/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/":$PYTHONPATH
export PYTHONPATH="/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/":$PYTHONPATH
export R_LIBS="/home/egusmao/R/x86_64-pc-linux-gnu-library/3.3:/projects/cheops.AE/software/R/R-3.3.3_intel_mkl/lib64/R/library":$R_LIBS

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/25_AB_Compartments/input_ab/"${SLURM_ARRAY_TASK_ID}".txt"
chromosome=`sed '1q;d' $inputFileName`
resolution=`sed '2q;d' $inputFileName`
chromSizesFileName=`sed '3q;d' $inputFileName`
enhChromSizesFileName=`sed '4q;d' $inputFileName`
eigenFileName=`sed '5q;d' $inputFileName`
dnaseFileName=`sed '6q;d' $inputFileName`
tempLocation=`sed '7q;d' $inputFileName`
outBedFileName=`sed '8q;d' $inputFileName`
outBwFileName=`sed '9q;d' $inputFileName`
outCountFileName=`sed '10q;d' $inputFileName`

# Execution
python /projects/ag-papan/eduardo/Papantonis_Stag/Code/25_AB_Compartments/2_formatAB.py $chromosome $resolution $chromSizesFileName $enhChromSizesFileName $eigenFileName $dnaseFileName $tempLocation $outBedFileName $outBwFileName $outCountFileName


