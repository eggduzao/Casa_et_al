
#Import
import os
import sys
from glob import glob

# Input
name = sys.argv[1]
pValue = sys.argv[2]
commandType = sys.argv[3]
treatmentFileName = sys.argv[4]
controlFileName = sys.argv[5]
tempLocation = sys.argv[6]
outputLocation = sys.argv[7]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

# Creating command
if(commandType == "Ed"):
  command = "macs2 callpeak -n "+name+" -t "+treatmentFileName+" -c "+controlFileName+" -p "+pValue+" -f BAM -g hs --keep-dup auto --outdir "+outputLocation+" --verbose 3 --nomodel --shift -100 --extsize 200 --nolambda --tempdir "+tempLocation+" --call-summits"
else:
  command = "macs2 callpeak -n "+name+" -t "+treatmentFileName+" -c "+controlFileName+" -p "+pValue+" -f BAM -g hs --outdir "+outputLocation+" --verbose 3 --tempdir "+tempLocation+" --call-summits"


# Applying MACS2
os.system(command)


