
rm(list = ls())

# Libraries
library(hicrep)

# Input
args = commandArgs(trailingOnly = TRUE)
matrixFileName1 = args[1]
matrixFileName2 = args[2]
corrFileName = args[3]
resolution = as.numeric(args[4])
smoothing = as.numeric(args[5])
maxDistInteract = as.numeric(args[6])

# Reading matrices
matrix1 = read.table(matrixFileName1, header = FALSE, sep = "\t")
matrix2 = read.table(matrixFileName2, header = FALSE, sep = "\t")

# Creating hicrep object
Pre_HiC <- prep(matrix1, matrix2, resolution, smoothing, maxDistInteract)

# Get the best smoothing parameter
#h_hat <- htrain(matrix1, matrix2, resolution, 5000000, 0:2)

# Performing test
SCC.out = get.scc(Pre_HiC, resolution, maxDistInteract)

# Writing corr (SCC score) / and SD (Standard deviation of SCC)
write(c(SCC.out[[3]], SCC.out[[4]]), file = corrFileName, ncolumns = 2, append = FALSE, sep = "\t")

rm(list = ls())
