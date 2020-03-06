
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

# Input
args <- commandArgs(trailingOnly = TRUE)
resolution = as.numeric(args[1])
logt = args[2]
domOrder = as.numeric(args[3])
maxDistInBin = as.numeric(args[4])
minD = as.numeric(args[5])
maxD = as.numeric(args[6])
minDP = as.numeric(args[7])
maxDP = as.numeric(args[8])
hthr = as.numeric(args[9])
t1thr = as.numeric(args[10])
inputFileName = args[11]
outputFilePrefix = args[12]

# Initialization
logt = FALSE
if(logt == "T"){logt = TRUE}

###################################################################################################
# Execution
###################################################################################################

# Reading input table
hic_mat = read.table(inputFileName, header = FALSE, sep = "\t")
hic_mat = hic_mat[,c(2,3,4)]
colnames(hic_mat) = c("n1", "n2", "count")

# Applying GMAP
res = rGMAP(hic_mat, resl = resolution, logt = logt, dom_order = domOrder, maxDistInBin = maxDistInBin, min_d = minD, max_d = maxD,
            min_dp = minDP, max_dp = maxDP, hthr = hthr, t1thr = t1thr)

# Writing results to CSV table
outputTADName = paste(outputFilePrefix,"_tad.txt",sep="")
write.table(res$tads, file = outputTADName, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)
outputHTADName = paste(outputFilePrefix,"_htad.txt",sep="")
write.table(res$hierTads, file = outputHTADName, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", row.names = FALSE, col.names = TRUE)


