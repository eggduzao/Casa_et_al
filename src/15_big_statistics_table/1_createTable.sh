#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=STAB
#SBATCH --output=STAB.out
#SBATCH --error=STAB.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=10:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=6-7

# Commands
# sbatch 1_createTable.sh
# squeue -u egusmao
# scancel 10159695

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
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/15_big_statistics_table/input_stab/"${SLURM_ARRAY_TASK_ID}".txt"
fileIsLoop=`sed '1q;d' $inputFileName`
stagFileName=`sed '2q;d' $inputFileName`
aliasFileName=`sed '3q;d' $inputFileName`
chromSizesFileName=`sed '4q;d' $inputFileName`
genomeFileName=`sed '5q;d' $inputFileName`
regionsFileName=`sed '6q;d' $inputFileName`
ensemblDictFileName=`sed '7q;d' $inputFileName`
enhancersFileName=`sed '8q;d' $inputFileName`
chrommHmmFileName=`sed '9q;d' $inputFileName`
expressionLabelList=`sed '10q;d' $inputFileName`
expressionFileNameList=`sed '11q;d' $inputFileName`
signalLabelList=`sed '12q;d' $inputFileName`
signalCountList=`sed '13q;d' $inputFileName`
signalFileNameList=`sed '14q;d' $inputFileName`
controlCountList=`sed '15q;d' $inputFileName`
controlFileNameList=`sed '16q;d' $inputFileName`
peakFileNameList=`sed '17q;d' $inputFileName`
motifLabelList=`sed '18q;d' $inputFileName`
motifFileNameList=`sed '19q;d' $inputFileName`
tempLocation=`sed '20q;d' $inputFileName`
outputFileName1=`sed '21q;d' $inputFileName`
outputFileName2=`sed '22q;d' $inputFileName`
outputFileName3=`sed '23q;d' $inputFileName`
outputFileName4=`sed '24q;d' $inputFileName`

# Execution
python /projects/ag-papan/eduardo/Papantonis_Stag/Code/15_big_statistics_table/1_createTable.py $fileIsLoop $stagFileName $aliasFileName $chromSizesFileName $genomeFileName $regionsFileName $ensemblDictFileName $enhancersFileName $chrommHmmFileName $expressionLabelList $expressionFileNameList $signalLabelList $signalCountList $signalFileNameList $controlCountList $controlFileNameList $peakFileNameList $motifLabelList $motifFileNameList $tempLocation $outputFileName1 $outputFileName2 $outputFileName3 $outputFileName4


