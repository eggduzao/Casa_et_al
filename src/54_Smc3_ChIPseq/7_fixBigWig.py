
#Import
import os
import sys
import numpy as np
from pysam import Samfile
from multiprocessing import Process

###################################################################################################
# Input
###################################################################################################

# Input
chromosome = sys.argv[1]
chromSizesFileName = sys.argv[2]
chromSizesFileEnhName = sys.argv[3]
mainBamFileName = sys.argv[4]
toAddBamFileNameList = sys.argv[5].split(",")
toRemoveBamFileNameList = sys.argv[6].split(",")
outWigFileName = sys.argv[7]

# Initialization
outLoc = "/".join(outWigFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outLoc
os.system(command)

###################################################################################################
# Functions
###################################################################################################

# Read chrom sizes
def get_chrom_sizes(genome_sizes_file_name):
  genome_sizes_dict = dict()
  genome_sizes_file = open(genome_sizes_file_name, "rU")
  for line in genome_sizes_file:
    ll = line.strip().split("\t")
    genome_sizes_dict[ll[0]] = int(ll[1])
  genome_sizes_file.close()
  chrom_list = sorted(genome_sizes_dict.keys())
  return chrom_list, genome_sizes_dict

# Compress vector
def compressVector(total_bins, vector):
  return [np.average(e) for e in np.array_split(np.array(vector), total_bins)]

# Fetch signal BAM wrapper
def wrapper_fetchSignalBam_function(args):
  fetchSignalBam(*args)

# Fetch signal BAM
def fetchSignalBam(total_bins, region, bamFile):
  vector = [0.0] * (region[2] - region[1])
  #try:
  for read in bamFile.fetch(region[0], region[1], region[2]):
    for k in range(max(read.reference_start, region[1]), min(read.reference_end, region[2])): vector[k - region[1]] += 1.0
  return compressVector(total_bins, vector)
  #except Exception: pass

def convert_to_bigwig(wig_file_name, genome_sizes_file_name, bw_file_name, remove_original = False):
  command = " ".join(["wigToBigWig", wig_file_name, genome_sizes_file_name, bw_file_name])
  os.system(command)
  if(remove_original): os.system("rm "+wig_file_name)

# Main
def fix_bigwig(chromosome, chromSizesFileName, chromSizesFileEnhName, mainBamFileName, toAddBamFileNameList, toRemoveBamFileNameList, outWigFileName):

  # Fixed parameters
  GENOME_WINDOW_SIZE = 1000000
  WINDOW_SIZE = 10
  TOTAL_BINS = GENOME_WINDOW_SIZE / WINDOW_SIZE

  # Get chrom sizes
  chrom_list, genome_sizes_dict = get_chrom_sizes(chromSizesFileName)
  
  # Open bam and wig files
  outWigFile = open(outWigFileName, "w")
  mainBamFile = Samfile(mainBamFileName, "rb")
  toAddBamFileList = [Samfile(e,"rb") for e in toAddBamFileNameList]
  toRemoveBamFileList = [Samfile(e,"rb") for e in toRemoveBamFileNameList]

  # Wig header
  wig_header = "fixedStep chrom="+chromosome+" start=1 step="+str(WINDOW_SIZE)
  outWigFile.write(wig_header+"\n")

  # Iterating on genomic regions for memory purposes
  for i in range(0, genome_sizes_dict[chromosome], GENOME_WINDOW_SIZE):

    # Region to fetch the signal
    region = [chromosome, i, min(i+GENOME_WINDOW_SIZE, genome_sizes_dict[chromosome])]

    # Fetch signals
    vector_list_add = []
    vector_list_rm = []
    mainSignal = fetchSignalBam(TOTAL_BINS, region, mainBamFile)
    for j in range(0, len(toAddBamFileList)):
      vector_list_add.append(fetchSignalBam(TOTAL_BINS, region, toAddBamFileList[j]))
    for j in range(0, len(toRemoveBamFileList)):
      vector_list_rm.append(fetchSignalBam(TOTAL_BINS, region, toRemoveBamFileList[j]))

    # Writing signals
    for j in range(0, TOTAL_BINS):
      vMain = mainSignal[j]
      vToAdd = sum([0.2 * e[j] for e in vector_list_add])
      vToRemove = sum([0.3 * e[j] for e in vector_list_rm])
      outWigFile.write(str(max((vMain + vToAdd) - vToRemove, 0.0))+"\n")

  # Termination
  mainBamFile.close()
  outWigFile.close()
  for e in toAddBamFileList: e.close()
  for e in toRemoveBamFileList: e.close()
  convert_to_bigwig(outWigFileName, chromSizesFileEnhName, ".".join(outWigFileName.split(".")[:-1] + ["bw"]), remove_original = False)

###################################################################################################
# Execution
###################################################################################################

# Fix bigwig
fix_bigwig(chromosome, chromSizesFileName, chromSizesFileEnhName, mainBamFileName, toAddBamFileNameList, toRemoveBamFileNameList, outWigFileName)


