
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys

# Input
maxDiffInBp = int(sys.argv[1])
tadFileName1 = sys.argv[2]
tadFileName2 = sys.argv[3]
tempLocation = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
chromosome = tadFileName1.split("/")[-2]
command = "mkdir -p "+tempLocation
os.system(command)
command = "mkdir -p "+"/".join(outputFileName.split("/")[:-1])+"/"
os.system(command)

###################################################################################################
# Functions
###################################################################################################

# Merge TADs with the same length
def mergeTad(tadFileName, tempLocation, newTadFileName):

  # Reading tad dictionary
  tadDict = dict()
  tadFile = open(tadFileName, "rU")
  tadFile.readline()
  for line in tadFile:
    ll = line.strip().split("\t")
    key = ":".join([ll[0], ll[1]])
    try:
      tad = tadDict[key]
      continue
    except Exception: tadDict[key] = [ll[0], ll[1]]
  tadFile.close()

  # Writing temporary file
  newTadFileTempName = tempLocation + "newTadFileTempName.txt"
  newTadFileTemp = open(newTadFileTempName, "w")
  for k in tadDict.keys():
    newTadFileTemp.write("\t".join(tadDict[k])+"\n")
  newTadFileTemp.close()

  # Sort file
  command = "sort -k1,1n "+newTadFileTempName+" > "+newTadFileName
  os.system(command)

# Each TAD = [[pos1, pos2], accountedFor]
def readTadList(tadFileName):
  tadList = []
  tadFile = open(tadFileName, "rU")
  for line in tadFile:
    ll = line.strip().split("\t")
    tadList.append([[int(ll[0]), int(ll[1])], False])
  tadFile.close()
  return tadList

def overlap(t1, t2):
  if t1[1] <= t2[0]: return -1 # interval1 is before interval2
  if t2[1] <= t1[0]: return 1 # interval1 is after interval2
  else: return 0 # interval1 overlaps interval2

###################################################################################################
# Read tad lists
###################################################################################################

# Merge equal TADs (equal in length)
newTadFileName1 = tempLocation + "newTadFileName1.txt"
newTadFileName2 = tempLocation + "newTadFileName2.txt"
mergeTad(tadFileName1, tempLocation, newTadFileName1)
mergeTad(tadFileName2, tempLocation, newTadFileName2)

# Reading TADs
tadList1 = readTadList(newTadFileName1)
tadList2 = readTadList(newTadFileName2)
len1 = len(tadList1)
len2 = len(tadList2)

# Creating dictionary of possible status
tadStatus = {"0": "equal", "1": "split", "2": "merge", "3": "shift_upstream", "4": "shift_downstream", "5": "new"}

###################################################################################################
# Create differential table
###################################################################################################

# Initialization
diffTable = [] # [[tad(s)1], [tad(s)2], status]
currTad1 = None
currTad2 = None
counter1 = 0
counter2 = 0

# Iterating through TADs
while True:

  # Fetching TAD1
  try:
    currTad1 = tadList1[counter1]
  # If there is no more TAD1
  except Exception:
    # If there are still TAD2
    if(counter2 < len2):
      # If TAD2 was already "accounted for" -> just move on
      if(currTad2[1]):
        counter2 += 1
      # If TAD2 was not "accounted for" -> insert as new and move on
      else:
        diffTable.append([["NA"], currTad2[0], 5])
        currTad2[1] = True
        counter2 += 1
      continue
    # If there are no more TAD2 -> BREAK
    else: break

  # Fetching TAD2
  try:
    currTad2 = tadList2[counter2]
  # If there is no more TAD2
  except Exception:
    # If there are still TAD1
    if(counter1 < len1):
      # If TAD1 was already "accounted for" -> just move on
      if(currTad1[1]):
        counter1 += 1
      # If TAD1 was not "accounted for" -> insert as new and move on
      else:
        diffTable.append([currTad1[0], ["NA"], 5])
        currTad1[1] = True
        counter1 += 1
      continue
    # If there are no more TAD1 -> BREAK
    else: break
  
  # Checking overlap status of the TADs
  ovl = overlap(currTad1[0], currTad2[0])

  # TAD1 is completely before TAD2
  if(ovl == -1):

    # If TAD1 was "accounted for" -> just move on
    if(currTad1[1]):
      counter1 += 1
      continue
    # If TAD1 was not "accounted for" -> insert as new and move on
    else:
      diffTable.append([currTad1[0], ["NA"], 5])
      currTad1[1] = True
      counter1 += 1

  # TAD1 is completely after TAD2
  elif(ovl == 1):

    # If TAD2 was "accounted for" -> just move on
    if(currTad2[1]):
      counter2 += 1
      continue
    # If TAD2 was not "accounted for" -> insert as new and move on
    else:
      diffTable.append([["NA"], currTad2[0], 5])
      currTad2[1] = True
      counter2 += 1

  # TADs overlap
  else:

    # Acceptance ranges
    aRange1S = [currTad1[0][0] - maxDiffInBp, currTad1[0][0] + maxDiffInBp]
    aRange1E = [currTad1[0][1] - maxDiffInBp, currTad1[0][1] + maxDiffInBp]
    aRange2S = [currTad2[0][0] - maxDiffInBp, currTad2[0][0] + maxDiffInBp]
    aRange2E = [currTad2[0][1] - maxDiffInBp, currTad2[0][1] + maxDiffInBp]

    # Variable contains whether the starting limits of TADs match
    sM = (aRange1S[0] <= currTad2[0][0] and aRange1S[1] >= currTad2[0][0]) or (aRange2S[0] <= currTad1[0][0] and aRange2S[1] >= currTad1[0][0])

    # Variable contains whether the ending limits of TADs match
    eM = (aRange1E[0] <= currTad2[0][1] and aRange1E[1] >= currTad2[0][1]) or (aRange2E[0] <= currTad1[0][1] and aRange2E[1] >= currTad1[0][1])

    # If both ends of the TAD match -> equal
    if(sM and eM):
      diffTable.append([currTad1[0], currTad2[0], 0])
      currTad1[1] = True
      currTad2[1] = True
      counter1 += 1
      counter2 += 1

    # If starting limits of TADs match
    elif(sM):
      
      # If TAD1 is longer -> split
      if(currTad1[0][1] > currTad2[0][1]):
        diffTable.append([currTad1[0], currTad2[0], 1])
        currTad1[1] = True
        currTad2[1] = True
        counter2 += 1

      # Else if TAD2 is longer -> merge
      elif(currTad1[0][1] < currTad2[0][1]):
        diffTable.append([currTad1[0], currTad2[0], 2])
        currTad1[1] = True
        currTad2[1] = True
        counter1 += 1

    # If ending limits of TADs match
    elif(eM):

      # If TAD1 is longer -> split
      if(currTad1[0][0] < currTad2[0][0]):
        diffTable.append([currTad1[0], currTad2[0], 1])
        currTad1[1] = True
        currTad2[1] = True
        counter1 += 1
        counter2 += 1

      # Else if TAD2 is longer -> merge
      elif(currTad1[0][0] > currTad2[0][0]):
        diffTable.append([currTad1[0], currTad2[0], 2])
        currTad1[1] = True
        currTad2[1] = True
        counter1 += 1
        counter2 += 1

    # If no limits of TADs match
    else:

      # If TAD1 is before TAD2 -> shift upstream
      if((currTad1[0][0] < currTad2[0][0]) and (currTad1[0][1] < currTad2[0][1])):
        diffTable.append([currTad1[0], currTad2[0], 3])
        currTad1[1] = True
        currTad2[1] = True
        counter1 += 1        

      # If TAD1 is after TAD2 -> shift downstream
      elif((currTad1[0][0] > currTad2[0][0]) and (currTad1[0][1] > currTad2[0][1])):
        diffTable.append([currTad1[0], currTad2[0], 4])
        currTad1[1] = True
        currTad2[1] = True
        counter2 += 1        

      # If TAD1 encompasses all TAD2 -> split
      elif((currTad1[0][0] < currTad2[0][0]) and (currTad1[0][1] > currTad2[0][1])):
        diffTable.append([currTad1[0], currTad2[0], 1])
        currTad1[1] = True
        currTad2[1] = True
        counter2 += 1        

      # If TAD1 is inside TAD2 -> merge
      elif((currTad1[0][0] > currTad2[0][0]) and (currTad1[0][1] < currTad2[0][1])):
        diffTable.append([currTad1[0], currTad2[0], 2])
        currTad1[1] = True
        currTad2[1] = True
        counter1 += 1        

###################################################################################################
# Writing differential table
###################################################################################################

# Iterating on diffTable and writing differential TADs
outputFile = open(outputFileName, "w")
for item in diffTable:
  try: tad1s = str(item[0][0])
  except Exception: tad1s = "NA"
  try: tad1e = str(item[0][1])
  except Exception: tad1e = "NA"
  try: tad2s = str(item[1][0])
  except Exception: tad2s = "NA"
  try: tad2e = str(item[1][1])
  except Exception: tad2e = "NA"
  try: stateNb = str(item[2])
  except Exception: stateNb = "NA"
  try: stateName = tadStatus[stateNb]
  except Exception: stateName = "NA"
  outputFile.write("\t".join([chromosome, tad1s, tad1e, tad2s, tad2e, stateName, stateNb])+"\n")
outputFile.close()

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)

# cut -f 1,2,3,6 /media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/3B9_5.txt | grep -v NA > ~/Desktop/3B9_5_1.bed
# cut -f 1,4,5,6 /media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/3B9_5.txt | grep -v NA > ~/Desktop/3B9_5_2.bed
# cut -f 1,2,3,6 /media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/69_127.txt | grep -v NA > ~/Desktop/69_127_1.bed
# cut -f 1,4,5,6 /media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/69_127.txt | grep -v NA > ~/Desktop/69_127_2.bed


