#!/bin/bash

# Slurm Parameters
#SBATCH --job-name=STAGME
#SBATCH --output=STAGME.out
#SBATCH --error=STAGME.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=30:00:00
#SBATCH --account=UniKoeln
#SBATCH --array=1-14

# Commands
# sbatch 1_performEnrichment.sh
# squeue -u egusmao
# scancel 9755292

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
export PYTHONPATH=$PYTHONPATH:"/home/egusmao/.local/lib/python2.7/site-packages"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/cutadapt-1.15/inst/lib/python2.7/site-packages/"
export PYTHONPATH=$PYTHONPATH:"/projects/ag-papan/install/MACS2-2.1.1.20160309/inst/lib/python2.7/site-packages/"

# Input
inputFileName="/projects/ag-papan/eduardo/Papantonis_Stag/Code/9_tf_enrichment/input/"${SLURM_ARRAY_TASK_ID}".txt"
fprThreshold=`sed '1q;d' $inputFileName`
pseudocounts=`sed '2q;d' $inputFileName`
multipleTestAlpha=`sed '3q;d' $inputFileName`
randProportion=`sed '4q;d' $inputFileName`
foregroundFileName=`sed '5q;d' $inputFileName`
backgroundFileName=`sed '6q;d' $inputFileName`
outputLocationMatching=`sed '7q;d' $inputFileName`
outputLocationEnrichment=`sed '8q;d' $inputFileName`
mkdir -p $outputLocationMatching
mkdir -p $outputLocationEnrichment

# If Regions vs Random analysis
if [ $randProportion = "." ]
then

  # Motif Matching - Foreground
  rgt-motifanalysis matching --organism hg19 --fpr $fprThreshold --pseudocounts $pseudocounts --output-location $outputLocationMatching --input-files $foregroundFileName

  # Motif Matching - Background
  rgt-motifanalysis matching --organism hg19 --fpr $fprThreshold --pseudocounts $pseudocounts --output-location $outputLocationMatching --input-files $backgroundFileName

  # Motif Enrichment
  rgt-motifanalysis enrichment --organism hg19 --matching-location $outputLocationMatching --multiple-test-alpha $multipleTestAlpha --output-location $outputLocationEnrichment --print-thresh "1.0" $backgroundFileName $foregroundFileName

else

  # Motif Matching - Foreground and Random
  rgt-motifanalysis matching --organism hg19 --fpr $fprThreshold --pseudocounts $pseudocounts --rand-proportion $randProportion --output-location $outputLocationMatching --input-files $foregroundFileName

  # Motif Enrichment
  rgt-motifanalysis enrichment --organism hg19 --matching-location $outputLocationMatching --multiple-test-alpha $multipleTestAlpha --output-location $outputLocationEnrichment --print-thresh "1.0" $outputLocationMatching"/random_regions.bed" $foregroundFileName

fi


