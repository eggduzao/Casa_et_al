
#Import
import os
import sys

# Input
bamFileNameList = sys.argv[1].split(",")
tempBamFileName = sys.argv[2]
outputBamFileName = sys.argv[3]

# Create output location
outputLoc = "/".join(outputBamFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLoc
os.system(command)

# Create temporary location
tempLoc = "/".join(tempBamFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+tempLoc
os.system(command)

# Execution
if(not bamFileNameList): print "Ta doido?"
if(len(bamFileNameList) == 1):
  command = "cp "+bamFileNameList[0]+" "+tempBamFileName
  os.system(command)
  command = "samtools sort -o "+outputBamFileName+" "+tempBamFileName
  os.system(command)
  command = "samtools index "+outputBamFileName
  os.system(command)
else:
  command = "samtools merge "+tempBamFileName+" "+" ".join(bamFileNameList)
  os.system(command)
  command = "samtools sort -o "+outputBamFileName+" "+tempBamFileName
  os.system(command)
  command = "samtools index "+outputBamFileName
  os.system(command)


