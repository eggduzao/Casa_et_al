
#Import
import os
import sys
from glob import glob

# Input
threshold = sys.argv[1]
dataType = sys.argv[2]
treatmentFileName = sys.argv[3]
controlFileName = sys.argv[4]
tempLocation = sys.argv[5]
outputLocation = sys.argv[6]

# Initialization
treatName = ".".join(treatmentFileName.split("/")[-1].split(".")[:-1])
tempFolder = tempLocation+treatName+"/"
command = "mkdir -p "+tempFolder
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

# Creating command
command = "macs2 callpeak -n "+treatName+" -t "+treatmentFileName+" -c "+controlFileName+" -p "+threshold+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --shift -100 --extsize 200 --nolambda --tempdir "+tempFolder+" --call-summits"

# Applying MACS2
os.system(command)


