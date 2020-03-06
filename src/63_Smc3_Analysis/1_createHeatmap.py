
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
ext = int(sys.argv[1])
zMax = sys.argv[2]
hheight = sys.argv[3]
featureSummitFileName = sys.argv[4]
bwFileName = sys.argv[5]
bwLabel = sys.argv[6]
tempLocation = sys.argv[7]
outputLocation = sys.argv[8]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Heatmaps
###################################################################################################

# Allowed chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Creating bed from peak files
featureSummitFile = open(featureSummitFileName,"r")
tempBedFileName = tempLocation+"bedfile.bed"
bedFile = open(tempBedFileName,"w")
for line in featureSummitFile:
  ll = line.strip().split("\t")
  if(ll[0] not in chrList): continue
  region = [ll[0], int(ll[1])-ext, int(ll[2])+ext]
  if(int(region[1]) < 0): continue
  bedFile.write("\t".join([str(e) for e in region])+"\n")
featureSummitFile.close()
bedFile.close()

# Creating heatmaps
tempMatFileName = tempLocation+"matrix.mat.gz"

# Creating matrix
command = "computeMatrix reference-point -S \""+bwFileName+"\" -R \""+tempBedFileName+"\" -a \""+str(ext)+"\" -b \""+str(ext)+"\" --referencePoint \"center\" --binSize \"1\" --sortRegions \"descend\" --missingDataAsZero --numberOfProcessors \"max/2\" --outFileName \""+tempMatFileName+"\""
os.system(command)

# Creating heatmap
command = "plotHeatmap -m \""+tempMatFileName+"\" -out \""+outputLocation+bwLabel+".png\" --heatmapHeight \""+hheight+"\" --zMax \""+zMax+"\" --dpi \"90\" --colorMap \"Oranges\" --missingDataColor \"white\" --refPointLabel \"Summit\" --yAxisLabel \""+bwLabel+" Signal\" --samplesLabel \""+bwLabel+"\" --legendLocation \"lower-center\" --plotFileFormat \"png\""
os.system(command)





