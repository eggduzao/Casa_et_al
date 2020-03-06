
###################################################################################################
# Import
###################################################################################################

# Library
rm(list=ls())
library(rGMAP)
set.seed(111)

# Parameters:
# hic_mat -      Either a 3 columns Hi-C contact matrix for a given chromosome, with each row
#                  corrsponding to the start bin, end bin and the contact number; or a n by n matrix,
#                  n is the number of bins for a given chromosome.
# resl -         The resolution (bin size), default 10kb.
# logt -         Do log-transformation or not, default TRUE.
# dom_order -    Maximum level of hierarchical structures, default 2 (call TADs and subTADs).
# maxDistInBin - Only consider contact whose distance is not greater than maxDistInBIn bins,
#                  default 200 bins (or 2Mb).
# min_d -        The minimum d (d: window size), default 25.
# max_d -        The maximum d (d: window size), default 100.
# min_dp -       The minmum dp (dp: lower bound of tad size), defalt 5.
# max_dp -       The maximum dp (dp: lower bound of tad size), defalt 10. min_d, max_d, and max_dp should be specified in number of bins.
# hthr -         The lower bound cutoff for posterior probability, default 0.95.
# t1thr -        Lower bound for t1 for calling TAD, default 0.5 quantile of test statistics of
#                TADs, 0.9 of subTADs.

# Value:
# A list includes following elements:
# tads -     A data frame with columns start, end indicates the start and end coordinates of each domain, respectively.
# hierTads - A data frame with columns start, end, dom_order, where dom_order indicates the hierarchical status of a domain,
#              1 indicates tads, 2 indicates subtads, and so on.
# params -   A data frame gives the final parameters for calling TADs

###################################################################################################
# Execution
###################################################################################################

# Input
ilf = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_matrix_files/10K_norm/"
olf = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/18_TADs/1_TADs_GMAP/"
inputList = c("3B9_5-", "3B9_5plus", "69_127-", "69_127plus")

# Input Loop
for(fname in inputList){

  # Chromosome
  chrList = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
              "chr21", "chr22", "chrX")

  # Chromosome Loop
  for(chrom in chrList){

    system(paste("echo",fname,chrom,sep = " "))

    # Parameters
    inputFileName = paste(ilf,fname,"/",chrom,"_bin.txt",sep="")
    outputFolder = paste(olf,fname,"/",sep="")
    system(paste("mkdir -p ",outputFolder,sep=""))
    outputTADName = paste(outputFolder,chrom,"_tad.txt",sep="")
    outputHTADName = paste(outputFolder,chrom,"_htad.txt",sep="")

    # Reading input table
    hic_mat = read.table(inputFileName, header = FALSE, sep = "\t")
    hic_mat = hic_mat[,c(2,3,4)]
    colnames(hic_mat) = c("n1", "n2", "count")

    # Applying GMAP
    resolution = 10000
    res = rGMAP(hic_mat, resl = resolution, logt = T, dom_order = 2, maxDistInBin = min(200, 2 * 10^6/resolution), min_d = 25, max_d = 100,
                min_dp = 5, max_dp = 10, hthr = 0.95, t1thr = 0.5)

    # Writing results to CSV table
    write.table(res$tads, file = outputTADName, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
    write.table(res$hierTads, file = outputHTADName, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)

  }

}


