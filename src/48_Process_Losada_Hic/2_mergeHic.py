
# Import
import os
import sys
from glob import glob

###################################################################################################
# INPUT
###################################################################################################

# Input
juicerCommandFile = sys.argv[1] # /projects/ag-papan/eduardo/juicer/scripts/common/mega.sh
genomeId = sys.argv[2] # hg19
restrictionEnzyme = sys.argv[3] # MboI
excludeFragments = sys.argv[4] # if 1 then -x (no FRAG is generated)
hicWorkFolderList = sys.argv[5].split(",")
tempLocation = sys.argv[6]
outputLocation = sys.argv[7] # ...../Results/10_Process_All_HiC_Data/2_Merged_Hic/XXXXX/

###################################################################################################
# EXECUTION
###################################################################################################

# Removing temp folder
command = "rm -rf "+tempLocation
os.system(command)

# Create temp and output location
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

# Create symbolic links to the tempLocation
for hicWorkFolder in hicWorkFolderList:
  command = "ln -s "+hicWorkFolder+" "+tempLocation
  os.system(command)

# Processing arguments
jcf = juicerCommandFile
gi = "-g "+genomeId
wf = "-d "+tempLocation
st = "-s "+restrictionEnzyme
if(restrictionEnzyme == "1"): nore = "-x"
else: nore = ""

# Creating command
command = " ".join([jcf, gi, wf, st, nore])
os.system(command)

# Postprocess
alignedFolder = tempLocation + "mega/aligned/"
for inFN in glob(alignedFolder+"*"):
  command = "mv "+inFN+" "+outputLocation
  os.system(command)

"""
[egusmao@cheops0 common]$ ./mega.sh -h
Usage: mega.sh -g genomeID [-d topDir] [-s site] [-h]
   genomeID must be defined in the script, e.g. "hg19" or "mm10" (default "hg19")
   [topDir] is the top level directory (default "/projects/ag-papan/eduardo/juicer/scripts/common") and must contain links to all merged_nodups files underneath it
   [site] must be defined in the script, e.g.  "HindIII" or "MboI" (default "MboI"); alternatively, this can be the restriction site file
* [stage]: must be one of "final", "postproc", or "early".\n    -Use "final" when the reads have been combined into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use "postproc" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use "early" for an early exit, before the final creation of the stats and\n     hic files
   -x: exclude fragment-delimited maps from Hi-C mega map (will run much faster)
   -h: print this help and exit
"""


