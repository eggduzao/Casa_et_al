
# Import
import os
import sys
from glob import glob

# Bam List
chromSizesFile = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.filter"
chromSizesFileEnh = "/usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced"
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/25_Smc3_ChIPseq/input/"
il1 = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/2_Merged_Bam_Files/"
il2 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/rao/bam/"
il3 = "/usr/users/egadegu/Projects/Wendt_Stag/Data/stag_bam_files/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/25_Smc3_ChIPseq/7_FixedBigWig/"
addBamFileNameList = [il2+"CHIP001_RAD21_untreated_SE.bam", il3+"STAG1_GFP.mrg-replicates.srt.bam", il3+"STAG2_GFP.mrg-replicates.srt.bam"]
remBamFileNameList = [il1+"Smc3_Input.bam", il2+"CHIP031_input_untreated_SE.bam"]

# Opening file
inputFileName = fl + "7_fbw.txt"
inputFile = open(inputFileName, "w")

# Chromosome List
chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Chromosome Loop
for chrom in chromList:

  # Input
  chromosome = chrom
  chromSizesFileName = chromSizesFile
  chromSizesFileEnhName = chromSizesFileEnh
  mainBamFileName = il1 + "Smc3_Untreated.bam"
  toAddBamFileNameList = ",".join(addBamFileNameList)
  toRemoveBamFileNameList = ",".join(remBamFileNameList)
  outWigFileName = ol + "Smc3_" + chrom + ".wig"

  # Execution
  inputFile.write(" ".join([chromosome, chromSizesFileName, chromSizesFileEnhName, mainBamFileName, toAddBamFileNameList, toRemoveBamFileNameList, outWigFileName])+"\n")

# Close
inputFile.close()

#cat Smc3_chr1.wig Smc3_chr2.wig Smc3_chr3.wig Smc3_chr4.wig Smc3_chr5.wig Smc3_chr6.wig Smc3_chr7.wig Smc3_chr8.wig Smc3_chr9.wig Smc3_chr10.wig Smc3_chr11.wig Smc3_chr12.wig Smc3_chr13.wig Smc3_chr14.wig Smc3_chr15.wig Smc3_chr16.wig Smc3_chr17.wig Smc3_chr18.wig Smc3_chr19.wig Smc3_chr20.wig Smc3_chr21.wig Smc3_chr22.wig Smc3_chrX.wig > Smc3.wig

#wigToBigWig Smc3.wig /usr/users/egadegu/rgtdata/hg19/chrom.sizes.hg19.enhanced Smc3.bw

