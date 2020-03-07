
# Import
import os
import sys

# Hic List
genomeID = "hg19"
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/39_Process_Rao_Matrix/input/"
il = "/usr/users/egadegu/Projects/Wendt_Stag/Data/rao_new_untreated_hic/"
tl = "/usr/users/egadegu/scratch/HTS/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/39_Process_Rao_Matrix/"
hicMat = [["Untreated"], ["Untreated_Synchronized"]]

# Opening file
inFileName = fl + "1_hts.txt"
inFile = open(inFileName,"w")

# Hic Loop
for hicList in hicMat:

  # Parameter
  outName = hicList[0]

  # Resolution List
  resList = ["5000", "10000", "25000", "50000", "100000", "250000", "500000"]

  # Resolution Loop
  for res in resList:

    # Input
    genome_id = genomeID
    resolution = res
    chrom_sizes_file_name = chromSizesFile
    input_hic_file_name_list = ",".join([il + e + ".hic" for e in hicList])
    temporary_location = tl + outName + "/" + res + "/"
    output_hic_file_name = ol + outName + "_" + res + ".hic"
    output_sparse_file_name = ol + outName + "_" + res + ".txt"

    # Creating files
    inFile.write(" ".join([genome_id, resolution, chrom_sizes_file_name, input_hic_file_name_list, temporary_location, output_hic_file_name, output_sparse_file_name])+"\n")

# Closing file
inFile.close()

