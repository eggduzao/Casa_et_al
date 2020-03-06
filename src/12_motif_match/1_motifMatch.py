
# Import
import os
import sys
from glob import glob

# Input
motifFileName = sys.argv[1]
regionsFileName = sys.argv[2]
tempLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Temporary Location
command = "mkdir -p "+tempLocation
os.system(command)

# Motif Matching
motifName = motifFileName.split("/")[-1].split(".")[0]
olMm = tempLocation+motifName+"MM/"
command = "rgt-motifanalysis matching --organism hg19 --fpr 0.001 --pseudocounts 1 --use-only-motifs "+motifFileName+" --output-location "+olMm+" --input-files "+regionsFileName
os.system(command)

# Moving to output file
outMMFileName = glob(olMm+"*_mpbs.bed")[0]
command = "mv "+outMMFileName+" "+outputFileName
os.system(command)

# Deleting files
command = "rm -rf "+tempLocation
os.system(command)


