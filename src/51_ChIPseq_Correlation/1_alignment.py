
#Import
import os
import sys
from glob import glob

# Input
analysisType = sys.argv[1]
minQuality = sys.argv[2]
ncores = sys.argv[3]
fastqGzFileName = sys.argv[4]
indexFileName = sys.argv[5]
tempLocation = sys.argv[6]
outputName = sys.argv[7]
outputLocation = sys.argv[8]

# Create temp and output location
command = "mkdir -p "+outputLocation
os.system(command)
tempFolder = tempLocation + outputName + "/"
command = "mkdir -p "+tempFolder
os.system(command)

# Paired end analysis
if(analysisType == "PE"):

  # Parameters
  fastqGzFileName = fastqGzFileName.split(",")

  # Fasta file(s)
  fastqFileName1 = tempFolder + outputName + "_R1.fa"
  fastqFileName2 = tempFolder + outputName + "_R2.fa"
  command = "gunzip -cd \""+fastqGzFileName[0]+"\" > "+fastqFileName1
  os.system(command)
  command = "gunzip -cd \""+fastqGzFileName[1]+"\" > "+fastqFileName2
  os.system(command)

  # Uncompressing index
  command = "unzip "+indexFileName+" -d "+tempFolder
  os.system(command)
  indexPrefixFileName = glob(tempFolder+"*.bt2")[0].split("/")[-1].split(".")[0]

  # Running alignment
  samFileName = tempFolder + outputName + ".sam"
  cwd = os.getcwd()
  command = "cd "+tempFolder+"; bowtie2 --no-mixed --no-discordant -X 1000 -p "+ncores+" -x "+indexPrefixFileName+" -1 "+fastqFileName1+" -2 "+fastqFileName2+" -S "+samFileName+"; cd "+cwd
  os.system(command)

  # Sam to Bam
  bamFileName = outputLocation + outputName + ".bam"
  command = "samtools view -bhS "+samFileName+" > "+bamFileName
  os.system(command)

  # Sorting, removing duplicates and indexing bam files
  on = ".".join(bamFileName.split(".")[:-1])
  qualityFileName = outputLocation + outputName + "_align_quality.txt"
  command = ("cd "+outputLocation+"; "
             "samtools sort "+on+".bam -o "+on+"_sorted.bam; "
             "samtools index "+on+"_sorted.bam; "
             "samtools rmdup -sS "+on+"_sorted.bam "+on+"_remdup.bam; "
             "samtools view -bq "+minQuality+" "+on+"_remdup.bam > "+on+".bam; "
             "samtools index "+on+".bam; "
             "rm "+on+"_sorted.bam "+on+"_remdup.bam "+on+"_sorted.bam.bai; "
             "samtools flagstat "+on+".bam > "+qualityFileName+" ")
  os.system(command)

# Single end analysis
if(analysisType == "SE"): 

  # Fasta file(s)
  fastqFileName = tempFolder + outputName + ".fa"
  command = "gunzip -cd \""+fastqGzFileName+"\" > "+fastqFileName
  os.system(command)

  # Uncompressing index
  command = "unzip "+indexFileName+" -d "+tempFolder
  os.system(command)
  indexPrefixFileName = glob(tempFolder+"*.bt2")[0].split("/")[-1].split(".")[0]

  # Perform trimming using trim galore
  command = "trim_galore --no_report_file --suppress_warn -q "+minQuality+" -o "+tempFolder+" "+fastqFileName
  os.system(command)
  trimmedFileName = glob(tempFolder+"*trimmed*")[0]

  # Running alignment
  samFileName = tempFolder + outputName + ".sam"
  cwd = os.getcwd()
  command = "cd "+tempFolder+"; bowtie2 --no-mixed --no-discordant -X 1000 -p "+ncores+" -x "+indexPrefixFileName+" -U "+trimmedFileName+" -S "+samFileName+"; cd "+cwd
  os.system(command)

  # Sam to Bam
  bamFileName = outputLocation + outputName + ".bam"
  command = "samtools view -bhS "+samFileName+" > "+bamFileName
  os.system(command)

  # Sorting, removing duplicates and indexing bam files
  on = ".".join(bamFileName.split(".")[:-1])
  qualityFileName = outputLocation + outputName + "_align_quality.txt"
  command = ("cd "+outputLocation+"; "
             "samtools sort "+on+".bam -o "+on+"_sorted.bam; "
             "samtools index "+on+"_sorted.bam; "
             "samtools rmdup -sS "+on+"_sorted.bam "+on+"_remdup.bam; "
             "samtools view -bq "+minQuality+" "+on+"_remdup.bam > "+on+".bam; "
             "samtools index "+on+".bam; "
             "rm "+on+"_sorted.bam "+on+"_remdup.bam "+on+"_sorted.bam.bai; "
             "samtools flagstat "+on+".bam > "+qualityFileName+" ")
  os.system(command)


