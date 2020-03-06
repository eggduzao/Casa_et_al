
# Import
import os
import sys

# Input
start = sys.argv[1]
end = sys.argv[2]
resolution = sys.argv[3]
name = sys.argv[4]
cond = sys.argv[5]
matrixFileName = sys.argv[6]
tadFileName = sys.argv[7]
newTadFileName = sys.argv[8]
outputFileName = sys.argv[9]

# Converting TADs
command = "python 2_format_gmap.py "+" ".join([start, end, tadFileName, newTadFileName])
os.system(command)

# Printing matrix
command = "R CMD BATCH '--args '"+resolution+"' '"+start+"' '"+name+"' '"+cond+"' '"+matrixFileName+"' '"+newTadFileName+"' '"+outputFileName+" 2_tadPlot.R 2_tadPlot.Rout"
os.system(command)

# Removing TEMP files
command = "rm "+" ".join([newTadFileName])
os.system(command)

#python 2_format_gmap.py 10000000 20000000 /media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/3B9_5-/chr1/T_2_80_25_100_5_10_0.95_0.5_htad.txt /home/egg/Projects/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/test.txt

#R CMD BATCH '--args 25000 10000000 name treat 400_800.mat test.txt res.pdf' 2_tadPlot.R 2_tadPlot.Rout

#R CMD BATCH '--args 25000 10000000 name treat 400_800.mat test.txt res.png' 2_tadPlot2.R 2_tadPlot2.Rout


