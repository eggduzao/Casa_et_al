
#Import
import os
import sys

# Input
sraFileName = sys.argv[1]
outputLocation = sys.argv[2]

# Execution
command = "mkdir -p "+outputLocation
os.system(command)
command = "fastq-dump.2.8.2 --split-3 --outdir "+outputLocation+" "+sraFileName
os.system(command)


