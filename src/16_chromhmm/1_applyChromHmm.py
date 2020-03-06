
# Import
import os
import sys

###################################################################################################
# Input
###################################################################################################

# Input
nStates = sys.argv[1]
nProc = sys.argv[2]
organism = sys.argv[3]
chromSizesFileName = sys.argv[4]
cellMarkFileName = sys.argv[5]
inputBamFolder = sys.argv[6]
tempLocation = sys.argv[7]
outputLocation = sys.argv[8]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Execution
###################################################################################################

# Creating binarized model
command = "java -mx16000M -jar /home/egg/Install/chromhmm-1.15/ChromHMM.jar BinarizeBam "+chromSizesFileName+" "+inputBamFolder+" "+cellMarkFileName+" "+tempLocation
os.system(command)

# Applying chromhmm
command = "java -mx16000M -jar /home/egg/Install/chromhmm-1.15/ChromHMM.jar LearnModel -p "+nProc+" "+tempLocation+" "+outputLocation+" "+nStates+" "+organism
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


