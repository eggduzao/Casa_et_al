
# Import
import os
import sys
from glob import glob

###################################################################################################
# INPUT
###################################################################################################

# Input
juicerCommandFile = sys.argv[1] # /projects/ag-papan/eduardo/juicer/scripts/juicer.sh
genomeId = sys.argv[2] # hg19
restrictionEnzyme = sys.argv[3] # MboI
expDescription = sys.argv[4] # -a 'This is a test experiment'
commandStage = sys.argv[5] # -S merge / dedup / final / postproc / . (for all)
useShort = sys.argv[6] # if "1" -r
chromSizesLocation = sys.argv[7] # -p /projects/ag-papan/eduardo/juicer/references/chrom.sizes.hg19
restrictionSiteFileName = sys.argv[8] # -y <FILE> ["." for none]
assemblyGenomeFileName = sys.argv[9] # -z /projects/ag-papan/eduardo/juicer/references/Homo_sapiens_assembly19.fa
juicerFolder = sys.argv[10] # -D /projects/ag-papan/eduardo/juicer/
numberOfThreads = sys.argv[11] # -t 8
fastq1FileName = sys.argv[12] 
fastq2FileName = sys.argv[13]
workingFolder = sys.argv[14] # -d /scratch/eduardo/JUICER/XXXXX/
outputLocation = sys.argv[15] # ...../Results/10_Process_All_HiC_Data/1_Juicer/XXXXX/

###################################################################################################
# EXECUTION
###################################################################################################

# Removing work folder
if(os.path.exists(workingFolder)):
  command = "rm -rf "+workingFolder
  os.system(command)

# Creating work folder
fastqFolder = workingFolder + "fastq/"
command = "mkdir -p "+fastqFolder
os.system(command)

# Create output location
command = "mkdir -p "+outputLocation
os.system(command)

# Copying fasta files to working folder
command = "cp "+fastq1FileName+" "+fastq2FileName+" "+fastqFolder
os.system(command)

# Converting names to contain _R1 and _R2
oldFile1Name = fastqFolder + fastq1FileName.split("/")[-1]
oldFile2Name = fastqFolder + fastq2FileName.split("/")[-1]
newFile1Name = fastqFolder + "HIC_R1.fastq.gz"
newFile2Name = fastqFolder + "HIC_R2.fastq.gz"
command = "mv "+oldFile1Name+" "+newFile1Name
os.system(command)
command = "mv "+oldFile2Name+" "+newFile2Name
os.system(command)

# Processing arguments
jcf = juicerCommandFile
gi = "-g "+genomeId
wf = "-d "+workingFolder
st = "-s "+restrictionEnzyme
ed = "-a '"+expDescription+"'"
if(commandStage == "."): cs = ""
else: cs = "-S "+commandStage
if(useShort == "1"): us = "-r"
else: us = ""
csl = "-p "+chromSizesLocation
if(restrictionEnzyme == "none"): rsf = ""
else: rsf = "-y "+restrictionSiteFileName
agf = "-z "+assemblyGenomeFileName
jf = "-D "+juicerFolder
noth = "-t "+numberOfThreads
if(restrictionEnzyme == "none"): nore = ""
else: nore = "-f"

# Creating command
command = " ".join([jcf, gi, wf, st, ed, cs, us, csl, rsf, agf, jf, noth, nore])
print command
os.system(command)

# Postprocess
alignedFolder = workingFolder + "aligned/"
for inFN in glob(alignedFolder+"*"):
  command = "mv "+inFN+" "+outputLocation
  os.system(command)

"""
Usage: juicer.sh [-g genomeID] [-d topDir] [-s site] [-a about] [-R end]
                 [-S stage] [-p chrom.sizes path] [-y restriction site file]
                 [-z reference genome file] [-D Juicer scripts directory]
                 [-b ligation] [-t threads] [-r] [-h] [-f] [-j]
* [genomeID] must be defined in the script, e.g. "hg19" or "mm10" (default 
  "hg19"); alternatively, it can be defined using the -z command
* [topDir] is the top level directory (default
  "/usr/users/egadegu/Projects/Papantonis_Temp/Code")
     [topDir]/fastq must contain the fastq files
     [topDir]/splits will be created to contain the temporary split files
     [topDir]/aligned will be created for the final alignment
* [site] must be defined in the script, e.g.  "HindIII" or "MboI" 
  (default "MboI")
* [about]: enter description of experiment, enclosed in single quotes
* -r: use the short read version of the aligner, bwa aln
  (default: long read, bwa mem)
* [end]: use the short read aligner on read end, must be one of 1 or 2 
* [chrom.sizes path]: enter path for chrom.sizes file
* [restriction site file]: enter path for restriction site file (locations of
  restriction sites in genome; can be generated with the script
  misc/generate_site_positions.py)
* [reference genome file]: enter path for reference sequence file, BWA index
  files must be in same directory
* [Juicer scripts directory]: set the Juicer directory,
  which should have scripts/ references/ and restriction_sites/ underneath it
  (default /usr/users/egadegu/Juicer/)
* [ligation junction]: use this string when counting ligation junctions
* [threads]: number of threads when running BWA alignment
* -f: include fragment-delimited maps in hic file creation
* -h: print this help and exit
"""


