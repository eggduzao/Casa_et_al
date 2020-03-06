
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from pysam import Samfile

# Input
resolution = int(sys.argv[1])
abTreatFileNameList = sys.argv[2].split(",")
abControlFileNameList = sys.argv[3].split(",")
treatMatrixFileName = sys.argv[4]
controlMatrixFileName = sys.argv[5]
treatTadFileName = sys.argv[6]
controlTadFileName = sys.argv[7]
temporaryLocation = sys.argv[8]
outputIntraFileName = sys.argv[9]
outputInterFileName = sys.argv[10]

# Initialization
command = "mkdir -p "+temporaryLocation
os.system(command)
outputLocation = "/".join(outputIntraFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def round_down(num, divisor):
  return num - (num%divisor)

def round_up(num, divisor):
  return num + (num%divisor)

def read_matrix_dictionary(matrix_file_name):

  # Initialization
  matrix_dict = dict()
  matrix_file = open(matrix_file_name,"rU")
  bad_list = ["NaN", "Inf", "-Inf", "NA", "0", "0.0"]

  # Populating matrix
  for line in matrix_file:
    ll = line.strip().split("\t")
    if(ll[3] in bad_list): continue
    p1 = int(ll[1])
    p2 = int(ll[2])
    minP = str(min(p1, p2))
    maxP = str(max(p1, p2))
    matrix_dict[":".join([ll[0], minP, maxP])] = float(ll[3])
    matrix_dict[":".join([ll[0], maxP, minP])] = float(ll[3])

  # Returning objects
  matrix_file.close()
  return matrix_dict

def summ_avg_tad_interaction(hic_dict, resolution, chrom, b1, b2):

  # Calculating sum
  summ = 0.0
  count = 0.0
  b1 = round_down(b1, 25000)
  b2 = round_up(b2, 25000)
  for i in range(b1, b2 - resolution, resolution):
    for j in range(i + resolution, b2, resolution):
      key = ":".join([chrom, str(i), str(j)])
      try:
        summ += hic_dict[key]
        count += 1.0
      except Exception: count += 1.0

  # Returning objects
  return summ, round(summ/count, 8)

def summ_intertad_interaction(hic_dict, resolution, chrom, t1_1, t1_2, t2_1, t2_2):

  # Calculating sum
  summ = 0.0
  count = 0.0
  t1_1 = round_down(t1_1, 25000)
  t1_2 = round_up(t1_2, 25000)
  t2_1 = round_down(t2_1, 25000)
  t2_2 = round_up(t2_2, 25000)
  for i in range(t1_1, t1_2, resolution):
    for j in range(t2_1, t2_2, resolution):
      key = ":".join([chrom, str(i), str(j)])
      try:
        summ += hic_dict[key]
        count += 1.0
      except Exception: count += 1.0

  # Returning objects
  return summ, round(summ/count, 8)

###################################################################################################
# Execution
###################################################################################################

# Read Matrices
treatMatrixDict = read_matrix_dictionary(treatMatrixFileName)
contrMatrixDict = read_matrix_dictionary(controlMatrixFileName)

# Open TADs
treatTadFile = Samfile(treatTadFileName, "rb")
controlTadFile = Samfile(controlTadFileName, "rb")

# Iterate through AB files
abList = ["AA", "AB", "BA", "BB"]
abIntraDict = dict([(e, []) for e in abList])
abInterDict = dict([(e, []) for e in abList])
for i in range(0, len(abTreatFileNameList)):

  # Initialization
  abTreatFileName = abTreatFileNameList[i]
  abControlFileName = abControlFileNameList[i]
  abTreatFile = open(abTreatFileName, "rU")
  abControlFile = open(abControlFileName, "rU")
  aaInter = 0.0; abInter = 0.0; baInter = 0.0; bbInter = 0.0
  aaIntra = 0.0; abIntra = 0.0; baIntra = 0.0; bbIntra = 0.0

  # Iterating on AB locations
  for tline in abTreatFile:

    # Initialization
    cline = abControlFile.readline()
    tt = tline.strip().split("\t")
    cc = cline.strip().split("\t")
    chrom = tt[0]; p1 = int(tt[1]); p2 = int(tt[2]); tcomp = tt[3]; ccomp = cc[3]

    # Fetching TADs
    tfetch = treatTadFile.fetch(chrom, p1, p2)
    ttadList = []
    for read in tfetch: ttadList.append([read.reference_start, read.reference_end])
    cfetch = treatTadFile.fetch(chrom, p1, p2)
    ctadList = []
    for read in cfetch: ctadList.append([read.reference_start, read.reference_end])

    # Iterating on TADs treatment
    treatIntraCount = 0.0
    treatInterCount = 0.0
    flagTtad = True
    try:
      prevTad = ttadList[0]
      prevSum, prevAvg = summ_avg_tad_interaction(treatMatrixDict, resolution, chrom, prevTad[0], prevTad[1])
      treatIntraCount += prevAvg
      for k in range(1, len(ttadList)):
        currTad = ttadList[k]
        currSum, currAvg = summ_avg_tad_interaction(treatMatrixDict, resolution, chrom, currTad[0], currTad[1])
        treatIntraCount += currAvg
        interSum, interAvg = summ_intertad_interaction(treatMatrixDict, resolution, chrom, prevTad[0], prevTad[1], currTad[0], currTad[1])
        treatInterCount += interAvg
        prevTad = [currTad[0], currTad[1]]
        prevSum, prevAvg = currSum, currAvg
    except Exception: flagTtad = False

    # Iterating on TADs control
    contrIntraCount = 0.0
    contrInterCount = 0.0
    flagCtad = True
    try:
      prevTad = ttadList[0]
      prevSum, prevAvg = summ_avg_tad_interaction(contrMatrixDict, resolution, chrom, prevTad[0], prevTad[1])
      contrIntraCount += prevAvg
      for k in range(1, len(ttadList)):
        currTad = ttadList[k]
        currSum, currAvg = summ_avg_tad_interaction(contrMatrixDict, resolution, chrom, currTad[0], currTad[1])
        contrIntraCount += currAvg
        interSum, interAvg = summ_intertad_interaction(contrMatrixDict, resolution, chrom, prevTad[0], prevTad[1], currTad[0], currTad[1])
        contrInterCount += interAvg
        prevTad = [currTad[0], currTad[1]]
        prevSum, prevAvg = currSum, currAvg
    except Exception: flagCtad = False

    if(not flagTtad or not flagCtad): continue

    if(ccomp == "A" and tcomp == "A"):
      try: aaIntra += ((treatIntraCount / len(ttadList)) / (contrIntraCount / len(ctadList)))
      except Exception: pass
      try: aaInter += ((treatInterCount / (len(ttadList)-1)) / (contrInterCount / (len(ctadList)-1)))
      except Exception: pass
    elif(ccomp == "A" and tcomp == "B"):
      try: abIntra += ((treatIntraCount / len(ttadList)) / (contrIntraCount / len(ctadList)))
      except Exception: pass
      try: abInter += ((treatInterCount / (len(ttadList)-1)) / (contrInterCount / (len(ctadList)-1)))
      except Exception: pass
    elif(ccomp == "B" and tcomp == "A"):
      try: baIntra += ((treatIntraCount / len(ttadList)) / (contrIntraCount / len(ctadList)))
      except Exception: pass
      try: baInter += ((treatInterCount / (len(ttadList)-1)) / (contrInterCount / (len(ctadList)-1)))
      except Exception: pass
    elif(ccomp == "B" and tcomp == "B"):
      try: bbIntra += ((treatIntraCount / len(ttadList)) / (contrIntraCount / len(ctadList)))
      except Exception: pass
      try: bbInter += ((treatInterCount / (len(ttadList)-1)) / (contrInterCount / (len(ctadList)-1)))
      except Exception: pass

  abIntraDict["AA"].append(aaIntra)
  abIntraDict["AB"].append(abIntra)
  abIntraDict["BA"].append(baIntra)
  abIntraDict["BB"].append(bbIntra)
  abInterDict["AA"].append(aaInter)
  abInterDict["AB"].append(abInter)
  abInterDict["BA"].append(baInter)
  abInterDict["BB"].append(bbInter)

treatTadFile.close()
controlTadFile.close()

# Writing output
maxV = max([max(len(abIntraDict[e]), len(abInterDict[e])) for e in abList])
outputIntraFile = open(outputIntraFileName, "w")
outputInterFile = open(outputInterFileName, "w")
outputIntraFile.write("\t".join(abList)+"\n")
outputInterFile.write("\t".join(abList)+"\n")
for i in range(0, maxV):
  vecIntra = []
  vecInter = []
  for ab in abList:
    try: vecIntra.append(str(abIntraDict[ab][i]))
    except Exception: vecIntra.append("NA")
    try: vecInter.append(str(abInterDict[ab][i]))
    except Exception: vecInter.append("NA")
  outputIntraFile.write("\t".join(vecIntra)+"\n")
  outputInterFile.write("\t".join(vecInter)+"\n")
outputIntraFile.close()
outputInterFile.close()


