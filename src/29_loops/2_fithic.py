
# Import
import os
import sys

# Input
resolution = sys.argv[1]
lowerBound = sys.argv[2]
upperBound = sys.argv[3]
equalOcBinNb = sys.argv[4]
splineTimes = sys.argv[5]
contactFileName = sys.argv[6]
fragmentFileName = sys.argv[7]
outputName = sys.argv[8]
outputLocation = sys.argv[9]

# Execution
command = "fithic -i "+contactFileName+" -f "+fragmentFileName+" -o "+outputLocation+" -r "+resolution+" -L "+lowerBound+" -U "+upperBound+" -p "+splineTimes+" -b "+equalOcBinNb+" -l "+outputName
os.system(command)

# Filtering output files
csFileName = outputLocation + outputName + ".fithic.log"
fitFileName = outputLocation + outputName + ".fithic_pass1.res25000.txt"
splineFileName = outputLocation + outputName + ".spline_pass1.res40000.significances.txt.gz"
newSplineFileName = outputLocation + outputName + ".txt.gz"

# Removing files
command = "rm "+" ".join([csFileName, fitFileName])
os.system(command)

# Renaming spline
command = "mv "+" ".join([splineFileName, newSplineFileName])
os.system(command)


