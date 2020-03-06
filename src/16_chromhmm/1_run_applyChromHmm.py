
# Import
import os
import sys

# Input
nStates = "20"
nProc = "4"
organism = "hg19"
chromSizesFileName = "/home/egg/Install/chromhmm-1.15/CHROMSIZES/hg19.txt"
cellMarkFileName = "/home/egg/Projects/Papantonis_Stag/Code/16_chromhmm/cellMarkFileTable.txt"
inputBamFolder = "/media/egg/New4TB/bamFolder/"
tempLocation = "/home/egg/Projects/Papantonis_Stag/Code/16_chromhmm/TEMP/"
outputLocation = "/home/egg/Projects/Papantonis_Stag/Results/16_chromhmm/"

# Execution
command = "python 1_applyChromHmm.py "+" ".join([nStates, nProc, organism, chromSizesFileName, cellMarkFileName, inputBamFolder, tempLocation, outputLocation])
os.system(command)


