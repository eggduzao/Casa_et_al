
###################################################################################################
# Import
###################################################################################################

# Library
rm(list=ls())
library(lattice)
library(reshape)
library(plotrix)
library(ggplot2)
library(MASS)
library(ggthemes)
library(gplots)
library(RColorBrewer)
library(plot3D)
set.seed(111)

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(signalName, listOfVectors, columnVec, outFileName){

  # Plot Parameters
  yMax = -1
  for(i in 1:length(listOfVectors)){
    listOfVectors[[i]] = listOfVectors[[i]][columnVec]
    value = max(listOfVectors[[i]])
    if(value > yMax){ yMax = value }
  }
  lenOfVec = length(listOfVectors[[1]])
  
  # Graph Parameters
  colVec = c("black", "firebrick3", "dodgerblue3")
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(1.0, length(listOfVectors))
  labVec = c("Intersection", "Plus Only", "Minus Only")

  # Initilize figure
  pdf(file = outFileName, width = 6, height = 5)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(0,yMax)
  plot(xRange, yRange, type="n", xlab="Distance from STAG1-Only Peak Center", ylab="Average Signal", main="", axes = FALSE, cex.lab = 1.0)

  # Axis
  #axis(side = 1, at = c(1,500,1000), labels = c("-5000", "Center", "+5000")) #columnVec = 1:1000
  #axis(side = 1, at = c(1,50,100), labels = c("-500", "Center", "+500")) # columnVec = 450:550
  #axis(side = 1, at = c(1,100,200), labels = c("-1000", "Center", "+1000")) # columnVec = 400:600
  axis(side = 1, at = c(1,150,300), labels = c("-1500", "Center", "+1500")) # columnVec = 350:650
  axis(side = 2, at = seq(from = 0.0, to = yMax, by = ceiling(yMax)/10))

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  # Legend
  par(xpd=TRUE)
  legend(lenOfVec*0.00, yMax+(yMax*0.25), title = paste(signalName, " Signal", sep=""), labVec, col = colVec, lty = ltyVec, lwd = lwdVec, ncol = 3, xpd = TRUE)

  # Termination
  dev.off()

}

# Binning
binFun <- function(vector){
  retVec = c()
  binSize = 10
  numberOfBins = floor(length(vector)/binSize)
  for(i in 1:numberOfBins){
    retVec = c(retVec, sum(vector[(((binSize*i)-binSize)+1):(((binSize*i)-binSize)+binSize)]))
  }
  return(retVec)
}

###################################################################################################
# Execution
###################################################################################################

#################################################
# Line Plot
#################################################

# Input
ifn = "/home/egg/Projects/Papantonis_Stag/Results/11_aggregate_plots/tables/"
ofn = "/home/egg/Projects/Papantonis_Stag/Results/11_aggregate_plots/graphs/"

datasetList = c("HCT116_DNase-seq_UW")

#datasetList = c("HCT116_ChIP-seq_H3K27ac_USC", "HCT116_ChIP-seq_H3K36me3_USC", "HCT116_ChIP-seq_H3K4me1_USC", "HCT116_ChIP-seq_H3K4me2_BROAD", "HCT116_ChIP-seq_H3K4me3_UW", "HCT116_ChIP-seq_H3K9ac_USC", "HCT116_ChIP-seq_H3K9me3_USC", "HCT116_ChIP-seq_H3K9me2_BROAD", "HCT116_ChIP-seq_H3K79me2_BROAD", "HCT116_ChIP-seq_H3K27me3_BROAD", "HCT116_ChIP-seq_H3K20me1_BROAD", "HCT116_ChIP-seq_H2AFZ_BROAD")

datasetList = c("HCT116_ChIP-seq_ATF3_HAIB", "HCT116_ChIP-seq_CBX3_HAIB", "HCT116_ChIP-seq_CEBPB_HAIB", "HCT116_ChIP-seq_CTCF_HAIB", "HCT116_ChIP-seq_CTCF_UW", "HCT116_ChIP-seq_EGR1_HAIB", "HCT116_ChIP-seq_ELF1_HAIB", "HCT116_ChIP-seq_EZH2_BROAD", "HCT116_ChIP-seq_EZH2phosphoT487_BROAD", "HCT116_ChIP-seq_FOSL1_HAIB", "HCT116_ChIP-seq_MAX_HAIB", "HCT116_ChIP-seq_POL2RA_USC", "HCT116_ChIP-seq_POLR2AphosphoS5_HAIB", "HCT116_ChIP-seq_REST_HAIB", "HCT116_ChIP-seq_SIN3A_HAIB", "HCT116_ChIP-seq_SP1_HAIB", "HCT116_ChIP-seq_TCF7L2_USC", "HCT116_ChIP-seq_TEAD4_HAIB", "HCT116_ChIP-seq_USF1_HAIB", "HCT116_ChIP-seq_YY1_HAIB", "HCT116_ChIP-seq_ZBTB33_HAIB", "HCT116_ChIP-seq_ZNF274_USC", "hg19_HCT116_ChIP-seq_SRF_HAIB", "HCT116_ChIP-seq_CTCF_BROAD")

inFileListI = c()
inFileListP = c()
inFileListM = c()
outNameList = c()
for(dataset in datasetList){
  inFileListI = c(inFileListI, paste(ifn,dataset,"__R_STAG1_RM_STAG1_intersection.bed.txt",sep=""))
  inFileListP = c(inFileListP, paste(ifn,dataset,"__R_STAG1_RM_STAG1_plusOnly.bed.txt",sep=""))
  inFileListM = c(inFileListM, paste(ifn,dataset,"__R_STAG1_RM_STAG1_minusOnly.bed.txt",sep=""))
  outNameList = c(outNameList, paste(ofn,dataset,".pdf",sep=""))
}

# Parameters
#columnVec = 1:1000
#columnVec = 450:550
#columnVec = 400:600
columnVec = 350:650

# Iterating on input
for(i in 1:length(inFileListI)){

  # Input structure
  listOfVectors = vector("list", 3)  

  # Reading table I
  inputTableI = read.table(inFileListI[i], header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inputTableI[,4:ncol(inputTableI)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectors[[1]] = signalVector  

  # Reading table P
  inputTableP = read.table(inFileListP[i], header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inputTableP[,4:ncol(inputTableP)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectors[[2]] = signalVector  

  # Reading table M
  inputTableM = read.table(inFileListM[i], header = FALSE, sep = "\t")
  signalVector = as.numeric(colMeans(inputTableM[,4:ncol(inputTableM)]))
  signalVector[is.na(signalVector)] = 0.0
  signalVector = binFun(signalVector)
  listOfVectors[[3]] = signalVector  

  # Line Plot
  linePlot(datasetList[i], listOfVectors, columnVec, outNameList[i])

}


