#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=GMAP
#SBATCH --output=GMAP.out
#SBATCH --error=GMAP.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=10:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-1472

# Commands
# sbatch 1_generateTAD_GMAP.sh
# squeue -u egusmao
# scancel 10010253

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
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/input/"${SLURM_ARRAY_TASK_ID}".txt"
resolution=`sed '1q;d' $inputFileName`
logt=`sed '2q;d' $inputFileName`
domOrder=`sed '3q;d' $inputFileName`
maxDistInBin=`sed '4q;d' $inputFileName`
minD=`sed '5q;d' $inputFileName`
maxD=`sed '6q;d' $inputFileName`
minDP=`sed '7q;d' $inputFileName`
maxDP=`sed '8q;d' $inputFileName`
hthr=`sed '9q;d' $inputFileName`
t1thr=`sed '10q;d' $inputFileName`
inFileName=`sed '11q;d' $inputFileName`
outputFilePrefix=`sed '12q;d' $inputFileName`

#tempDir="/scratch/eduardo/1/"
#mkdir -p $tempDir
#cp /projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/1_generateTAD_GMAP.R $tempDir
#cd $tempDir
#R CMD BATCH '--args 25000 T 2 80 25 100 5 10 0.95 0.5 /projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/3B9_5-/chr1_bin.txt /projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/res' 1_generateTAD_GMAP.R 1_generateTAD_GMAP.Rout
#mv $tempDir"res"* /projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/

# Execution
tempDir="/scratch/eduardo/"${SLURM_ARRAY_TASK_ID}"/"
mkdir -p $tempDir
cp /projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/1_generateTAD_GMAP.R $tempDir
cd $tempDir
R CMD BATCH '--args '$resolution' '$logt' '$domOrder' '$maxDistInBin' '$minD' '$maxD' '$minDP' '$maxDP' '$hthr' '$t1thr' '$inFileName' '$outputFilePrefix 1_generateTAD_GMAP.R 1_generateTAD_GMAP.Rout


