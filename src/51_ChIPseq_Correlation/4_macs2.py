
#Import
import os
import sys
from glob import glob

# Input
dataType = sys.argv[1]
treatmentFileName = sys.argv[2]
controlFileName = sys.argv[3]
tempLocation = sys.argv[4]
outputLocation = sys.argv[5]

# Initialization
treatName = ".".join(treatmentFileName.split("/")[-1].split(".")[:-1])
tempFolder = tempLocation+treatName+"/"
command = "mkdir -p "+tempFolder
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

# Creating command
if(dataType == "DNASE"):
  if(controlFileName.split("/")[-1].split(".")[0] == "N"): command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --nolambda --tempdir "+tempFolder+" --call-summits"
  else: command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -c "+controlFileName+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --nolambda --tempdir "+tempFolder+" --call-summits"
elif(dataType == "CHIP"):
  if(controlFileName.split("/")[-1].split(".")[0] == "N"): command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --shift -100 --extsize 200 --nolambda --tempdir "+tempFolder+" --call-summits"
  else: command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -c "+controlFileName+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --shift -100 --extsize 200 --nolambda --tempdir "+tempFolder+" --call-summits"
elif(dataType == "PE"):
  if(controlFileName.split("/")[-1].split(".")[0] == "N"): command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -f BAMPE -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --tempdir "+tempFolder+" --call-summits"
  else: command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -c "+controlFileName+" -f BAMPE -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --tempdir "+tempFolder+" --call-summits"
else: print "choose a data type"

# Applying MACS2
os.system(command)


