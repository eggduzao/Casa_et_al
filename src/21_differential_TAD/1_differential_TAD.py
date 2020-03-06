
# Import
import os
import sys

# Input
resolution = int(sys.argv[1])
maxDiffBin = int(sys.argv[2])
tadFileName1 = sys.argv[3]
tadFileName2 = sys.argv[4]
tempLocation = sys.argv[5]
outputFileName1 = sys.argv[6]
outputFileName2 = sys.argv[7]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+"/".join(outputFileName1.split("/")[:-1])+"/"
os.system(command)

###################################################################################################
# Functions
###################################################################################################

# Create tad list
def tadList(inFileName):
  retList = []
  inFile = open(inFileName, "rU")
  header = "\t".join(inFile.readline().strip().split("\t")[:-1])
  for line in inFile:
    ll = line.strip().split("\t")
    retList.append(":".join(ll[:-1]))
  inFile.close()
  return header, retList

# Create expanded tad list
def expTadList(inFileName, resolution, maxDiffBin):
  retList = []
  inFile = open(inFileName, "rU")
  header = "\t".join(inFile.readline().strip().split("\t")[:-1])
  for line in inFile:
    ll = line.strip().split("\t")
    p1 = int(ll[0]); p2 = int(ll[1])
    for i in range(p1-(maxDiffBin*resolution), p1+((maxDiffBin+1)*resolution), resolution):
      for j in range(p2-(maxDiffBin*resolution), p2+((maxDiffBin+1)*resolution), resolution):
        retList.append(":".join([str(i),str(j)]))
  inFile.close()
  return header, retList

###################################################################################################
# Get original, extended and differential lists
###################################################################################################

# Get lists
header, originalList1 = tadList(tadFileName1)
header, originalList2 = tadList(tadFileName2)
header, extendedList1 = expTadList(tadFileName1, resolution, maxDiffBin)
header, extendedList2 = expTadList(tadFileName2, resolution, maxDiffBin)

diffList1 = list(set(originalList1) - set(extendedList2)) 
diffList2 = list(set(originalList2) - set(extendedList1))

#outFile = open("a.txt","w")
#outFile.write("\n".join(originalList1))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > originalList1.txt")

#outFile = open("a.txt","w")
#outFile.write("\n".join(originalList2))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > originalList2.txt")

#outFile = open("a.txt","w")
#outFile.write("\n".join(extendedList1))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > extendedList1.txt")

#outFile = open("a.txt","w")
#outFile.write("\n".join(extendedList2))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > extendedList2.txt")

#outFile = open("a.txt","w")
#outFile.write("\n".join(diffList1))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > diffList1.txt")

#outFile = open("a.txt","w")
#outFile.write("\n".join(diffList2))
#outFile.close()
#os.system("sort -t\":\" -k1,1n a.txt  > diffList2.txt")

#os.system("rm a.txt")

###################################################################################################
# Writing differential TADs
###################################################################################################

# Write list 1
tempFileName = tempLocation+"tempFileName.bed"
tempFile = open(tempFileName, "w")
for diff in diffList1: tempFile.write("\t".join(diff.split(":"))+"\n")
tempFile.close()

# Sort list 1
command = "sort -k3,3n -k1,1n "+tempFileName+" > "+outputFileName1
os.system(command)

# Inserting header 1
command = "sed -i -e '1i"+header+"\\' "+outputFileName1
os.system(command)

#####

# Write list 2
tempFileName = tempLocation+"tempFileName.bed"
tempFile = open(tempFileName, "w")
for diff in diffList2: tempFile.write("\t".join(diff.split(":"))+"\n")
tempFile.close()

# Sort list 2
command = "sort -k3,3n -k1,1n "+tempFileName+" > "+outputFileName2
os.system(command)

# Inserting header 2
command = "sed -i -e '1i"+header+"\\' "+outputFileName2
os.system(command)

#####

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


