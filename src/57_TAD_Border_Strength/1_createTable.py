
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
from glob import glob

# Input
resolution = int(sys.argv[1])
matrixFileName = sys.argv[2]
tadFileName = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
outputLocation = "/".join(outputFileName.split("/")[:-1])
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

# Reading hic matrix
hic_dict = read_matrix_dictionary(matrixFileName)

# Reading first TAD
tad_file = open(tadFileName, "rU")
ll = tad_file.readline().strip().split("\t")
prevTad = [ll[0], int(ll[1]), int(ll[2])]
prevSumm, prevTadAvg = summ_avg_tad_interaction(hic_dict, resolution, prevTad[0], prevTad[1], prevTad[2])
nextTad = None; nextTadSumm = 0.0; nextTadAvg = 0.0
interTadSumm = 0.0; interTadAvg = 0.0

# Iterating on TAD files and writing distribution
outputFile = open(outputFileName, "w")
for line in tad_file:
  ll = line.strip().split("\t")
  if(ll[0] != prevTad[0]):
    prevTad = [ll[0], int(ll[1]), int(ll[2])]
    continue
  nextTad = [ll[0], int(ll[1]), int(ll[2])]
  nextTadSumm, nextTadAvg = summ_avg_tad_interaction(hic_dict, resolution, nextTad[0], nextTad[1], nextTad[2])
  interTadSumm, interTadAvg = summ_intertad_interaction(hic_dict, resolution, nextTad[0], prevTad[1], prevTad[2], nextTad[1], nextTad[2])
  outputFile.write(str(((prevTadAvg+nextTadAvg+1)/2)/(interTadAvg+1))+"\n")
  prevTad = [ll[0], int(ll[1]), int(ll[2])]
  prevSumm = nextTadSumm; prevTadAvg = nextTadAvg
tad_file.close()
outputFile.close()


