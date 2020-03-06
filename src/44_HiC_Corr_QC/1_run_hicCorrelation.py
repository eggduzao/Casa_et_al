
# Import
import os
import sys

# Fixed Parameters
juicerCommand = "juicertools"
kindOfMatrix = "observed"
kindOfNormalization = "NONE"
unitOfResolution = "BP"
resolution = "25000"
smoothing = "2"
maxDistInteract = "1000000"
regionsFileName = "/home/egg/Projects/hic_corr/data/splitRegionsAll.txt"

# Variable parameters
p1 = "/home/egg/Projects/hic_corr/data/yulia_data/"
s1 = "_30.hic"
p2 = "/home/egg/Projects/hic_corr/data/rao_data/GSE63525_"
s2 = "_combined_30.hic"
hicFileList = [p1+"75726"+s1, p1+"75728"+s1, p1+"75729"+s1, p1+"79643"+s1, p1+"79644"+s1,
               p1+"79645"+s1, p1+"79646"+s1, p2+"IMR90"+s2]

# Execution
for i in range(0,len(hicFileList)-1):
  for j in range(i+1,len(hicFileList)):
    hicFileName1 = hicFileList[i]
    hicFileName2 = hicFileList[j]
    name1 = hicFileName1.split("/")[-1].split("_")[0]
    name2 = hicFileName2.split("/")[-1].split("_")[0]
    if("IMR90" in hicFileName1): name1 = "IMR90"
    if("IMR90" in hicFileName2): name2 = "IMR90"
    expName = "_".join([name1, name2])
    outputFileName="/home/egg/Projects/hic_corr/result_all/"+expName+".txt"
    command = "./run.sh "+" ".join([expName, juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution,
              resolution, smoothing, maxDistInteract, regionsFileName, hicFileName1, hicFileName2, outputFileName])
    os.system(command)


