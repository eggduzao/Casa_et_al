
#Import
import os
import sys
from glob import glob

# Input
alignType = sys.argv[1] # SE or PE
minQuality = sys.argv[2].strip()
ncores = sys.argv[3].strip()
if(alignType == "SE"): fastqGzFileName = sys.argv[4].strip()
elif(alignType == "PE"): fastqGzFileName = sys.argv[4].strip().split(",")
indexFileName = sys.argv[5].strip()
tempLocation = sys.argv[6].strip()
outputLocation = sys.argv[7].strip()

if(alignType == "SE"):

  # Initialization
  #fastqName = ".".join(fastqGzFileName.split("/")[-1].split(".")[:-2])
  #tempFolder = tempLocation+fastqName+"/"
  #command = "mkdir -p "+tempFolder
  #os.system(command)
  #fastqFileName = tempFolder+".".join(fastqGzFileName.split("/")[-1].split(".")[:-1])
  #command = "gunzip -c \""+fastqGzFileName+"\" > "+fastqFileName
  #os.system(command)

  # Initialization
  fastqName = ".".join(fastqGzFileName.split("/")[-1].split(".")[:-1])
  tempFolder = tempLocation+fastqName+"/"
  command = "mkdir -p "+tempFolder
  os.system(command)
  fastqFileName = fastqGzFileName

  # Uncompressing index
  command = "unzip "+indexFileName+" -d "+tempFolder
  os.system(command)
  indexPrefixFileName = glob(tempFolder+"*.bt2")[0].split("/")[-1].split(".")[0]

  # Perform trimming using trim galore
  command = "trim_galore --no_report_file --suppress_warn -q "+minQuality+" -o "+tempFolder+" "+fastqFileName
  os.system(command)
  trimmedFileName = glob(tempFolder+"*trimmed*")[0]

  # Running alignment
  ff = ".".join(fastqGzFileName.split("/")[-1].split(".")[:-1]).split("_")
  fastqName = "_".join(ff[:5])
  samFileName = tempFolder+fastqName+".sam"
  cwd = os.getcwd()
  command = "cd "+tempFolder+"; bowtie2 -X 2000 -p "+ncores+" -x "+indexPrefixFileName+" -U "+trimmedFileName+" -S "+samFileName+"; cd "+cwd
  os.system(command)

  # Sam to Bam
  bamFileName = outputLocation+fastqName+".bam"
  command = "samtools view -bhS "+samFileName+" > "+bamFileName
  os.system(command)

  # Sorting, removing duplicates and indexing bam files
  on = ".".join(bamFileName.split(".")[:-1])
  qualityFileName = outputLocation+fastqName+"_align_quality.txt"
  command = ("cd "+outputLocation+"; "
             "samtools sort "+on+".bam -o "+on+"_sorted.bam; "
             "samtools index "+on+"_sorted.bam; "
             "samtools rmdup -sS "+on+"_sorted.bam "+on+"_remdup.bam; "
             "samtools view -bq "+minQuality+" "+on+"_remdup.bam > "+on+".bam; "
             "samtools index "+on+".bam; "
             "rm "+on+"_sorted.bam "+on+"_remdup.bam "+on+"_sorted.bam.bai; "
             "samtools flagstat "+on+".bam > "+qualityFileName+" ")
  os.system(command)

if(alignType == "PE"):

  # Initialization
  #fastqName1 = ".".join(fastqGzFileName[0].split("/")[-1].split(".")[:-2])
  #fastqName2 = ".".join(fastqGzFileName[1].split("/")[-1].split(".")[:-2])
  #tempFolder = tempLocation+fastqName1+"/"
  #command = "mkdir -p "+tempFolder
  #os.system(command)
  #fastqFileName1 = tempFolder+".".join(fastqGzFileName[0].split("/")[-1].split(".")[:-1])
  #fastqFileName2 = tempFolder+".".join(fastqGzFileName[1].split("/")[-1].split(".")[:-1])
  #command = "gunzip -c \""+fastqGzFileName[0]+"\" > "+fastqFileName1
  #os.system(command)
  #command = "gunzip -c \""+fastqGzFileName[1]+"\" > "+fastqFileName2
  #os.system(command)

  # Initialization
  fastqName1 = ".".join(fastqGzFileName[0].split("/")[-1].split(".")[:-1])
  fastqName2 = ".".join(fastqGzFileName[1].split("/")[-1].split(".")[:-1])
  tempFolder = tempLocation+fastqName1+"/"
  command = "mkdir -p "+tempFolder
  os.system(command)
  fastqFileName1 = fastqGzFileName[0]
  fastqFileName2 = fastqGzFileName[1]

  # Uncompressing index
  command = "unzip "+indexFileName+" -d "+tempFolder
  os.system(command)
  indexPrefixFileName = glob(tempFolder+"*.bt2")[0].split("/")[-1].split(".")[0]

  # Perform trimming using trim galore
  #command = "trim_galore --no_report_file --suppress_warn -q "+minQuality+" -o "+tempFolder+" "+fastqFileName1
  #os.system(command)
  #trimmedFileName1 = glob(tempFolder+"*_PE_*_1_*trimmed*")[0]
  #command = "trim_galore --no_report_file --suppress_warn -q "+minQuality+" -o "+tempFolder+" "+fastqFileName2
  #os.system(command)
  #trimmedFileName2 = glob(tempFolder+"*_PE_*_2_*trimmed*")[0]

  # Running alignment
  ff = ".".join(fastqGzFileName[0].split("/")[-1].split(".")[:-1]).split("_")
  fastqName = "_".join(ff[:5])
  samFileName = tempFolder+fastqName+".sam"
  cwd = os.getcwd()
  command = "cd "+tempFolder+"; bowtie2 --no-mixed --no-discordant -X 1000 -p "+ncores+" -x "+indexPrefixFileName+" -1 "+fastqFileName1+" -2 "+fastqFileName2+" -S "+samFileName+"; cd "+cwd
  os.system(command)

  # Sam to Bam
  bamFileName = outputLocation+fastqName+".bam"
  command = "samtools view -bhS "+samFileName+" > "+bamFileName
  os.system(command)

  # Sorting, removing duplicates and indexing bam files
  on = ".".join(bamFileName.split(".")[:-1])
  qualityFileName = outputLocation+fastqName+"_align_quality.txt"
  command = ("cd "+outputLocation+"; "
             "samtools sort "+on+".bam -o "+on+"_sorted.bam; "
             "samtools index "+on+"_sorted.bam; "
             "samtools rmdup -sS "+on+"_sorted.bam "+on+"_remdup.bam; "
             "samtools view -bq "+minQuality+" "+on+"_remdup.bam > "+on+".bam; "
             "samtools index "+on+".bam; "
             "rm "+on+"_sorted.bam "+on+"_remdup.bam "+on+"_sorted.bam.bai; "
             "samtools flagstat "+on+".bam > "+qualityFileName+" ")
  os.system(command)


