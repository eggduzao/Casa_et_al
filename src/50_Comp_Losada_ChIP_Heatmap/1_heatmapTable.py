
# Import
import os
import sys
import math
import pyBigWig
from pysam import Samfile
from numpy import mean, array_split

###################################################################################################
# INPUT
###################################################################################################

# Input
numberOfBins = int(sys.argv[1])
totalExtension = int(sys.argv[2])
totalReadsW = int(sys.argv[3])
totalReadsL = int(sys.argv[4])
bedFileNameW = sys.argv[5]
bedFileNameL = sys.argv[6]
bamFileNameW = sys.argv[7]
bamFileNameL = sys.argv[8]
tempLocation = sys.argv[9]
outputFileName = sys.argv[10]

###################################################################################################
# FUNCTIONS
###################################################################################################

def file_len(fname):
  i = -1
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
  return i + 1

def fetch_counts(downstream_extension, upstream_extension, value_to_add, input_file, region, filetype = "bam"):

  # Creating vector
  vector = []

  if(filetype == "bam"):

    # Creating vector
    total_length = (region[2]-region[1])
    incr = 2 * (downstream_extension + upstream_extension)
    vector = [0.0] * (total_length+(2*incr))
    region1i = region[1] - incr
    region2i = region[2] + incr
    if(region1i < 0): return None

    # Fetching bam signal
    for read in input_file.fetch(region[0], region1i, region2i):
      if(read.is_reverse):
        read_start = max(read.reference_end - downstream_extension, region1i)
        read_end = min(read.reference_end + upstream_extension, region2i)
      else:
        read_start = max(read.reference_start - upstream_extension, region1i)
        read_end = min(read.reference_start + downstream_extension, region2i)
      for i in range(read_start, read_end): vector[i-region1i] += value_to_add

    # Trimming vector
    vector = vector[(incr):(len(vector)-2*incr)]
    for i in range(0, len(vector)):
      if(not vector or math.isnan(vector[i]) or math.isinf(vector[i])): vector[i] = 0.0

  elif(filetype == "bw"):

    # Creating vector
    total_length = (region[2]-region[1])
    vector = [0.0] * total_length

    # Fetching bw signal
    try: valuesVec = input_file.values(region[0], region[1], region[2])
    except Exception: return vector
    for i in range(0, len(valuesVec)):
      if(not valuesVec[i] or math.isnan(valuesVec[i]) or math.isinf(valuesVec[i])): vector[i] = 0.0
      else: vector[i] = valuesVec[i]

  # Returning objects
  return vector

def fetch_binarized_vector(number_of_bins, signal_vector):

  # Returning objetcs
  try:
    return [mean(e) if sum(e) > 0 else 0.0 for e in array_split(signal_vector, number_of_bins)]
  except Exception:
    return None

def cut_file(inputFile, outputFile):

  command = "cut -f 1,2,3 "+inputFile+" > "+outputFile
  os.system(command)

def sort_file(inputFile, outputFile):

  command = "sort -k1,1 -k2,2n "+inputFile+" > "+outputFile
  os.system(command)

def extend_from_summit(extension, inputFile, outputFile):

  inFile = open(inputFile, "rU")
  outFile = open(outputFile, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid = (int(ll[1]) + int(ll[2])) / 2
    p1 = mid - extension
    p2 = mid + extension
    outFile.write("\t".join([ll[0], str(p1), str(p2)])+"\n")
  inFile.close()
  outFile.close()

def average_peaks(extension, inputFile, outputFile):

  inFile = open(inputFile, "rU")
  outFile = open(outputFile, "w")
  for line in inFile:
    ll = line.strip().split("\t")
    mid1 = (int(ll[1]) + int(ll[2])) / 2
    mid2 = (int(ll[4]) + int(ll[5])) / 2
    p1 = mid1 - extension
    p2 = mid1 + extension
    p3 = mid2 - extension
    p4 = mid2 + extension
    outFile.write("\t".join([ll[0], str(p1), str(p2), ll[0], str(p3), str(p4)])+"\n")
  inFile.close()
  outFile.close()

def create_files(inputFile1Name, inputFile2Name, temporaryLocation, outputFileName):

  cutFileName1 = temporaryLocation + "cutFileName1.bed"
  cut_file(inputFile1Name, cutFileName1)

  ext200FileName1 = temporaryLocation + "ext100FileName1.bed"
  extend_from_summit(100, cutFileName1, ext200FileName1)

  sortFileName1 = temporaryLocation + "sortFileName1.bed"
  sort_file(ext200FileName1, sortFileName1)

  cutFileName2 = temporaryLocation + "cutFileName2.bed"
  cut_file(inputFile2Name, cutFileName2)

  ext200FileName2 = temporaryLocation + "ext200FileName2.bed"
  extend_from_summit(100, cutFileName2, ext200FileName2)

  sortFileName2 = temporaryLocation + "sortFileName2.bed"
  sort_file(ext200FileName2, sortFileName2)

  tempIntFileName = temporaryLocation + "tempIntFileName.bed"
  command = "intersectBed -wa -wb -a "+sortFileName1+" -b "+sortFileName2+" > "+tempIntFileName
  os.system(command)

  avgIntFileName = temporaryLocation + "avgIntFileName.bed"
  average_peaks(100, tempIntFileName, avgIntFileName)

  tempOnly1FileName = temporaryLocation + "tempOnly1FileName.bed"
  command = "intersectBed -v -wa -a "+sortFileName1+" -b "+sortFileName2+" > "+tempOnly1FileName
  os.system(command)

  tempOnly2FileName = temporaryLocation + "tempOnly2FileName.bed"
  command = "intersectBed -v -wa -a "+sortFileName2+" -b "+sortFileName1+" > "+tempOnly2FileName
  os.system(command)

  #catFileName = temporaryLocation + "catFileName.bed"
  command = "cat " + " ".join([avgIntFileName, tempOnly1FileName, tempOnly2FileName]) + " > " + outputFileName
  os.system(command)

  #sort_file(catFileName, outputFileName)

  return file_len(avgIntFileName), file_len(tempOnly1FileName), file_len(tempOnly2FileName)

###################################################################################################
# EXECUTION
###################################################################################################

# Create locations
command = "mkdir -p "+tempLocation
os.system(command)
outputLocation = "/".join(outputFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outputLocation
os.system(command)

# Scaling factor
rpmW = 1.0 / (totalReadsW / 1000000.0)
rpmL = 1.0 / (totalReadsL / 1000000.0)

# Fetching peak file
mergedPeakFileName = tempLocation + "mergedPeakFileName.bed"
nInt, n1Only, n2Only = create_files(bedFileNameW, bedFileNameL, tempLocation, mergedPeakFileName)

# Writing signal to table
bamFileW = Samfile(bamFileNameW, "rb")
bamFileL = Samfile(bamFileNameL, "rb")
inFile = open(mergedPeakFileName, "rU")
outputFile = open(outputFileName, "w")
counter = 0
for line in inFile:
  ll = line.strip().split("\t")
  chrom = ll[0]
  if(len(ll) > 4):
    p1 = int(ll[1]) - (totalExtension/2); p2 = int(ll[2]) + (totalExtension/2); p3 = int(ll[4]) - (totalExtension/2); p4 = int(ll[4]) + (totalExtension/2)
  else:
    p1 = int(ll[1]) - (totalExtension/2); p2 = int(ll[2]) + (totalExtension/2); p3 = int(ll[1]) - (totalExtension/2); p4 = int(ll[2]) + (totalExtension/2)
  signal1 = fetch_binarized_vector(numberOfBins, fetch_counts(100, 100, rpmW, bamFileW, [chrom, p1, p2], "bam"))
  signal2 = fetch_binarized_vector(numberOfBins, fetch_counts(100, 100, rpmL, bamFileL, [chrom, p3, p4], "bam"))
  if(counter < nInt): rowname = "1"
  elif(counter < (nInt+n1Only)): rowname = "2"
  else: rowname = "3"
  counter += 1
  if(signal1 and signal2): outputFile.write("\t".join([rowname] + [str(e) for e in signal1 + signal2])+"\n")
bamFileW.close()
bamFileL.close()
inFile.close()
outputFile.close()


