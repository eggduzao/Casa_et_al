
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
import math
import numpy as np

# Input
juicerCommand = sys.argv[1]
kindOfMatrix = sys.argv[2]
kindOfNormalization = sys.argv[3]
unitOfResolution = sys.argv[4]
resolution = sys.argv[5]
chromSizesFileName = sys.argv[6]
inputHicFileName = sys.argv[7]
tempLocation = sys.argv[8]
outputMatrixFileName = sys.argv[9]
outputScoreFileName = sys.argv[10]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)
outputLocation = "/".join(outputMatrixFileName.split("/")[:-1])
command = "mkdir -p "+outputLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def hic_to_sparse(juicer_command, kind_of_matrix, kind_of_normalization, unit_of_resolution, resolution, chrom_sizes_file_name, input_hic_file_name, temp_location, output_matrix_file_name, output_score_file_name):

  # Chromosome sizes dict
  chromSizeDict = dict()
  chromSizesFile = open(chrom_sizes_file_name,"rU")
  for line in chromSizesFile:
    ll = line.strip().split("\t")
    chromSizeDict[ll[0]] = ll[1]
  chromSizesFile.close()

  # Opening output file to merge everything on the go
  outputMatrixFile = open(output_matrix_file_name, "w")
  scoreVector = []

  # Output location for individual chromosomes
  outputLocation = output_matrix_file_name.split(".")[0] + "/"
  command = "mkdir -p "+outputLocation
  os.system(command)

  # Creating individual sparse chromosome matrices
  for chrom in sorted(chromSizeDict.keys()):

    # Initialization
    chrWoChr = chrom.split("chr")[-1]
    region = ":".join([chrWoChr,"1",str(chromSizeDict[chrom])])
    outputChromFileName = outputLocation + chrom + ".txt"
    outputChromFile = open(outputChromFileName, "w")

    # Creating sparse matrix
    tempFileName = temp_location + chrom + ".txt"
    command = " ".join([juicer_command, "dump", kind_of_matrix, kind_of_normalization, input_hic_file_name, region, region, unit_of_resolution, resolution, tempFileName])
    os.system(command)

    # Writing entries
    tempFile = open(tempFileName, "rU")
    for line in tempFile:
      ll = line.strip().split("\t")
      value = float(ll[2])
      if(math.isnan(value) or not np.isfinite(value)): continue
      scoreVector.append(value)
      outputMatrixFile.write("\t".join([chrom]+ll)+"\n")
      outputChromFile.write("\t".join([chrom]+ll)+"\n")
    tempFile.close()
    outputChromFile.close()

  # Calculating percentiles
  outputScoreFile = open(output_score_file_name, "w")
  scoreVectorNp = np.array(scoreVector)
  for i in range(0,101):
    outputScoreFile.write("\t".join([str(i), str(np.percentile(scoreVectorNp, i))])+"\n")
  outputScoreFile.close()

  # Termination
  outputMatrixFile.close()

###################################################################################################
# Execution
###################################################################################################

hic_to_sparse(juicerCommand, kindOfMatrix, kindOfNormalization, unitOfResolution, resolution, chromSizesFileName, inputHicFileName, tempLocation, outputMatrixFileName, outputScoreFileName)


