
# Import
import os
import sys
from glob import glob

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
fitFileName = glob(outputLocation + outputName + ".fithic_pass1.res*.txt")[0]
splineFileName = glob(outputLocation + outputName + ".spline_pass1.res*.significances.txt.gz")[0]
newSplineFileName = outputLocation + outputName + ".txt.gz"
#finalSplineFileName = outputLocation + outputName + ".txt"

# Removing files
command = "rm "+" ".join([csFileName, fitFileName])
os.system(command)

# Renaming spline
command = "mv "+" ".join([splineFileName, newSplineFileName])
os.system(command)

# Unziping file
#command = "gzip -cd "+newSplineFileName+" > "+finalSplineFileName
#os.system(command)

# Removing gz file
#command = "rm "+newSplineFileName
#os.system(command)


