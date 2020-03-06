
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math

# Input
loopBins = int(sys.argv[1])
resolution = int(sys.argv[2])
regionsFileName = sys.argv[3]
hicFile1Name = sys.argv[4]
hicFile2Name = sys.argv[5]
outputFileName = sys.argv[6]

# Parameters
promExt = 2000

###################################################################################################
# Functions
###################################################################################################

def read_regions(regions_file_name):

  # Reading regions for loop-o-gram
  regions = []
  inputFile = open(regions_file_name,"rU")
  for line in inputFile:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = float(ll[1]); p2 = float(ll[2])
    r1 = str(int(math.floor(p1 / float(resolution)) * float(resolution)) - ((loopBins/2) * int(resolution)))
    r2 = str(int(math.ceil((p1+1) / float(resolution)) * float(resolution)) + ((loopBins/2) * int(resolution)))
    r3 = str(int(math.floor(p2 / float(resolution)) * float(resolution)) - ((loopBins/2) * int(resolution)))
    r4 = str(int(math.ceil((p2+1) / float(resolution)) * float(resolution)) + ((loopBins/2) * int(resolution)))
    region = [chrom, r1, r2, chrom, r3, r4]
    regions.append(region)
  inputFile.close()

  # Returning objects
  return regions

def read_hic_dictionary(hic_file_name):

  # Reading HiC matrix
  hic_dict = dict()
  matrixFile = open(hic_file_name,"rU")
  for line in matrixFile:
    ll = line.strip().split("\t")
    hic_dict[":".join(ll[:3])] = float(ll[3])
  matrixFile.close()

  # Returning objects
  return hic_dict

def update_loop_o_gram(loop_o_gram, region1, region2, hic_dict1, hic_dict2, resolution):

  # Updating loop-o-gram
  counterI = 0
  chrom = region1[0]
  for i in range(int(region1[1]),int(region1[2]),resolution):
    counterJ = 0
    for j in range(int(region2[1]),int(region2[2]),resolution):
      try:
        value1 = float(hic_dict1[":".join([chrom,str(i),str(j)])])
        value2 = float(hic_dict2[":".join([chrom,str(i),str(j)])])
        loop_o_gram[counterI][counterJ] = loop_o_gram[counterI][counterJ] + (value1/value2)
      except Exception: pass
      counterJ += 1
    counterI += 1
    
def creating_loop_o_gram_one_region(loop_bins, resolution, regions_vec, hic_dict1, hic_dict2):

  # Iterating on regions
  counter = 0.0
  loop_o_gram = [[0.0] * (loop_bins+1) for e in range(0,(loop_bins+1))]
  for region in regions_vec:

    # Regions
    region1 = region[:3]
    region2 = region[3:]

    # Updating matrix
    update_loop_o_gram(loop_o_gram, region1, region2, hic_dict1, hic_dict2, resolution)
 
    # Update average counter
    counter += 1.0

  # Making average
  for i in range(0,len(loop_o_gram)):
    for j in range(0,len(loop_o_gram[i])):
      loop_o_gram[i][j] = loop_o_gram[i][j] / counter

  # Returning objects
  return loop_o_gram

def writing_loop_o_gram(loop_o_gram, output_file_name):

  # Writing loop-o-gram to file
  outputMatrixFile = open(output_file_name, "w")
  for vec in loop_o_gram:
    toPrint = []
    for v in vec: toPrint.append(str(v))
    outputMatrixFile.write("\t".join(toPrint)+"\n")
  outputMatrixFile.close()

def create_table(loop_bins, resolution, regions_file_name, hic_file_1_name, hic_file_2_name, output_file_name):

  # Read regions
  regions_vec = read_regions(regions_file_name)

  # Read HiC signal dictionary
  hic_dict1 = read_hic_dictionary(hic_file_1_name)
  hic_dict2 = read_hic_dictionary(hic_file_2_name)

  # Creating loop-o-gram
  loop_o_gram = creating_loop_o_gram_one_region(loop_bins, resolution, regions_vec, hic_dict1, hic_dict2)

  # Writing loop-o-gram
  writing_loop_o_gram(loop_o_gram, output_file_name)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_table(loopBins, resolution, regionsFileName, hicFile1Name, hicFile2Name, outputFileName)


