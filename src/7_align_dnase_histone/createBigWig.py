
#Import
import os
import sys
from ireg.alignment import BamFile

###################################################################################################
# Input
###################################################################################################

# Input
downstreamExt = int(sys.argv[1])
upstreamExt = int(sys.argv[2])
genomeSizesFileName = sys.argv[3]
bamFileName = sys.argv[4]
bwFileName = sys.argv[5]

###################################################################################################
# Execution
###################################################################################################

# Parameters
forwardShift = 0
reverseShift = 0

# Creating bigwig
bamFile = BamFile(bamFileName, treat_as = "file")
bamFile.create_bigwig(downstreamExt, upstreamExt, forwardShift, reverseShift, genomeSizesFileName, bwFileName, remove_wig = True)
bamFile.close()


