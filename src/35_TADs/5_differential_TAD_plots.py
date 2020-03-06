
# Import
import os
import sys

# Input
start = sys.argv[1]
end = sys.argv[2]
resolution = sys.argv[3]
name1 = sys.argv[4]
cond1 = sys.argv[5]
name2 = sys.argv[6]
cond2 = sys.argv[7]
matrixFileName1 = sys.argv[8]
matrixFileName2 = sys.argv[9]
tadFileName1 = sys.argv[10]
tadFileName2 = sys.argv[11]
diffTadFileName1 = sys.argv[12]
diffTadFileName2 = sys.argv[13]
newTadFileName1 = sys.argv[14]
newTadFileName2 = sys.argv[15]
newDiffTadFileName1 = sys.argv[16]
newDiffTadFileName2 = sys.argv[17]
outputFileName = sys.argv[18]

# Converting TADs
command = "python 2_format_gmap.py "+" ".join([start, end, tadFileName1, newTadFileName1])
os.system(command)
command = "python 2_format_gmap.py "+" ".join([start, end, tadFileName2, newTadFileName2])
os.system(command)
command = "python 2_format_gmap.py "+" ".join([start, end, diffTadFileName1, newDiffTadFileName1])
os.system(command)
command = "python 2_format_gmap.py "+" ".join([start, end, diffTadFileName2, newDiffTadFileName2])
os.system(command)

# Printing matrix
command = "R CMD BATCH '--args '"+resolution+"' '"+start+"' '"+name1+"' '"+cond1+"' '"+name2+"' '"+cond2+"' '"+matrixFileName1+"' '"+matrixFileName2+"' '"+newTadFileName1+"' '"+newTadFileName2+"' '"+newDiffTadFileName1+"' '"+newDiffTadFileName2+"' '"+outputFileName+" 2_differential_TAD_plots.R 2_differential_TAD_plots.Rout"
os.system(command)

# Merging TADs
command = "montage "+outputFileName+"_1.pdf "+outputFileName+"_2.pdf -tile 2x1 -geometry +0+0 "+outputFileName+".pdf"
#command = "pdfjam "+outputFileName+"_1.pdf "+outputFileName+"_2.pdf --nup 2x1 --landscape --outfile "+outputFileName+".pdf"
os.system(command)

# Removing TEMP files
command = "rm "+" ".join([newTadFileName1, newTadFileName2, newDiffTadFileName1, newDiffTadFileName2, outputFileName+"_1.pdf", outputFileName+"_2.pdf"])
os.system(command)


