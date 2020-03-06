
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
plusRegionFileName = sys.argv[1]
minusRegionFileName = sys.argv[2]
tempLocation = sys.argv[3]
plusOnlyOutputFileName = sys.argv[4]
minusOnlyOutputFileName = sys.argv[5]
intersectionOutputFileName = sys.argv[6]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Preparing files
###################################################################################################

# Remove score column
tempCutPlusFileName = tempLocation+"cutplus.txt"
command = "cut -f 1,2,3,4,5,6 "+plusRegionFileName+" > "+tempCutPlusFileName
os.system(command)

tempCutMinusFileName = tempLocation+"cutminus.txt"
command = "cut -f 1,2,3,4,5,6 "+minusRegionFileName+" > "+tempCutMinusFileName
os.system(command)

# Sorting files
tempSortPlusFileName = tempLocation+"sortplus.txt"
command = "sort "+tempCutPlusFileName+" > "+tempSortPlusFileName
os.system(command)

tempSortMinusFileName = tempLocation+"sortminus.txt"
command = "sort "+tempCutMinusFileName+" > "+tempSortMinusFileName
os.system(command)

###################################################################################################
# Intersections
###################################################################################################

# Only plus
tempOnlyPlusFileName = tempLocation+"temponlyplus.txt"
command = "comm -23 "+tempSortPlusFileName+" "+tempSortMinusFileName+" > "+tempOnlyPlusFileName
os.system(command)
command = "sort -k1,1 -k2,2n "+tempOnlyPlusFileName+" > "+plusOnlyOutputFileName
os.system(command)

# Only minus
tempOnlyMinusFileName = tempLocation+"temponlyminus.txt"
command = "comm -13 "+tempSortPlusFileName+" "+tempSortMinusFileName+" > "+tempOnlyMinusFileName
os.system(command)
command = "sort -k1,1 -k2,2n "+tempOnlyMinusFileName+" > "+minusOnlyOutputFileName
os.system(command)

# Intersection
tempIntFileName = tempLocation+"tempint.txt"
command = "comm -12 "+tempSortPlusFileName+" "+tempSortMinusFileName+" > "+tempIntFileName
os.system(command)
command = "sort -k1,1 -k2,2n "+tempIntFileName+" > "+intersectionOutputFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


