#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=LGP
#SBATCH --output=LGP.%A_%a.out
#SBATCH --error=LGP.%A_%a.err
#SBATCH --nodes=1
#SBATCH --mem=60gb
#SBATCH --time=72:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=56-60

# Commands
# sbatch 2_metaGenePlot.sh
# squeue -u egusmao
# scancel 11012821

# Modules
module add python/2.7.5-V2
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
export R_LIBS="/home/egusmao/R/x86_64-pc-linux-gnu-library/3.3:/projects/cheops.AE/software/R/R-3.3.3_intel_mkl/lib64/R/library":$R_LIBS

# Input
inputFileName="/projects/ag-papan/eduardo/Wendt_Stag/Code/10_Losada_Meta_Plots_Genes/input/2_lgp.txt"
parameters=`sed "${SLURM_ARRAY_TASK_ID}q;d" $inputFileName`

# Execution
tempDir="/scratch/egusmao/LGP/"${SLURM_ARRAY_TASK_ID}"/"
mkdir -p $tempDir
cp /projects/ag-papan/eduardo/Wendt_Stag/Code/10_Losada_Meta_Plots_Genes/2_metaGenePlot.R $tempDir
cd $tempDir
Rscript --verbose 2_metaGenePlot.R $parameters > 2_metaGenePlot.Rout


