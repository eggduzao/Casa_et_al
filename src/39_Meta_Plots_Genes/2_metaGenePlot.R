
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
library(OneR)
set.seed(111)

# Input
args <- commandArgs(trailingOnly = TRUE)
label = as.character(args[1])
inputTableFileName = as.character(args[2])
outputAggrFileName = as.character(args[3])
outputHeatFileName = as.character(args[4])
outputClusFileName = as.character(args[5])

###################################################################################################
# Functions
###################################################################################################

# Line plot
linePlot <- function(yMax, yMin, legendList, listOfVectors, numberOfElements, outFileName){
  
  # Graph Parameters
  lenOfVec = length(listOfVectors[[1]])
  colVec = c("black")
  ltyVec = rep(1, length(listOfVectors))
  lwdVec = rep(3.0, length(listOfVectors))
  labVec = legendList

  # Initilize figure
  pdf(file = outFileName, width = 10, height = 8)
  par(mar = c(5,5,4,2))

  # Initialize plot
  xRange = c(1,lenOfVec)
  yRange = c(yMin,yMax)
  plot(xRange, yRange, type="n", xlab="Genomic Position", ylab="Average Signal Count (RPM)", main=labVec, axes = FALSE, cex.main = 2.0, cex.lab = 1.5)

  # Axis
  axis(side = 1, at = c(1,500,1000,1500,2000,2500), labels = c("TSS-6000", "TSS-3000", "TSS", "TTS", "TTS+3000", "TTS+6000"), cex.axis=1.2)
  axis(side = 2, at = round(seq(from = yMin, to = yMax, length.out = 5), 2), cex.axis=1.2)

  # Lines
  xVecLines = seq(from = 1, to = lenOfVec)
  for(i in 1:length(listOfVectors)){
    lines(xVecLines, listOfVectors[[i]], type = "l", lwd = lwdVec[i], lty = ltyVec[i], col = colVec[i])
  }

  grid()
  abline(v = 1000, lty = 2, col = "darkred")
  abline(v = 1500, lty = 2, col = "darkred")

  text(25, y = yMin + (0.95*(yMax-yMin)), paste("n = ",numberOfElements,sep=""), font = 3, cex = 1.2)

  # Termination
  dev.off()

}

# Heatmap
createHeatmap <- function(table, plotLabel, nClusters, outFileName, outClusterFileName, printRowNames = FALSE){

  # Parameters
  graphWidth = 8
  graphHeight = 10
  myDist = function(p1) dist(p1, method="euclidean")
  myHclust = function(p2) hclust(p2, method="ward.D")
  heatmapMargins = c(1,3)
  heatmapLmat =  rbind(c(5,5,4,0), c(3,2,1,0))
  heatmapLwid = c(1.5, 6.0, 0.5, 2.0)
  heatmapLhei = c(0.2, 4)
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "row"
  heatmapRowv = TRUE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapCexRow = 1.5
  heatmapOffsetRow = -0.3
  heatmapLabCol = FALSE
  heatmapSepColor = "black"
  heatmapColSep = c(0, 100, 200, 300, 400, 500)
  heatmapSepWidth = c(0.01, 0.01)
  heatmapRowSep = c(0, nrow(table))

  # Minimum and maximum values
  minV = as.numeric(quantile(table, c(.05)))
  maxV = as.numeric(quantile(table, c(.95)))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(minV, maxV,length=101))
  # Perform clustering and save in a file
  etable = dist(table, method = "euclidean")
  fit = hclust(etable, method="ward.D")
  groups = cutree(fit, k=nClusters)
  groupNames  = names(groups)
  groups = as.numeric(groups)
  clusterColors = brewer.pal(12, "Set3")
  clusterVec = clusterColors[groups]

  # Changing table row names
  spl = strsplit(rownames(table), ":", fixed = TRUE)
  newRowNames = c()
  for(i in 1:nrow(table)){
    newRowNames = c(newRowNames, spl[[i]][4])
  }
  rownames(table) = newRowNames
  if(printRowNames == TRUE){
    heatmapLabRow = rownames(table)
  } else{
    heatmapLabRow = FALSE
  }

  # Initializing graph
  png(filename = outFileName, width = 600, height = 800)
  par()

  # Heatmap
  out = heatmap.2(table, col = colorScheme, breaks = hmbreaks, distfun = myDist, hclustfun = myHclust, margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize, trace = heatmapTrace,
            tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv,
            key=heatmapKey, cexRow = heatmapCexRow, offsetRow = heatmapOffsetRow, labRow = heatmapLabRow, labCol = heatmapLabCol, colsep = heatmapColSep, rowsep = heatmapRowSep, lmat = heatmapLmat,
            RowSideColors = clusterVec)

  # Heatmap colorkey
  colkey(col = colorScheme, clim = c(minV, maxV), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL, side = 4, length = 0.70, width = 1.0, dist = -0.1, shift = -0.2, line = NA, pos = NA,
         outer = FALSE, font = NA, lty = 1, lwd = 1, lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 0.3, padj = NA,
         cex.axis = 1.0, mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Cluster colorkey
  colkey(col = clusterColors, clim = c(1, 13), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL, side = 4, length = 0.4, width = 1.0, dist = -0.1, shift = 0.28, addlines = FALSE,
         breaks = NULL, at = 1:13, labels = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", ""), tick = FALSE, line = NA,
         pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 2.2,
         padj = -0.5, cex.axis = 1.0, mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Heatmap annotation
  par(xpd = "NA")
  text(0.38, 1.07, plotLabel, adj = NULL, cex = 1.5)
  text(0.06, 1.045, "TSS-6000", adj = NULL, cex = 0.8)
  text(0.19, 1.045, "TSS-3000", adj = NULL, cex = 0.8)
  text(0.32, 1.045, "TSS", adj = NULL, cex = 0.8)
  text(0.45, 1.045, "TTS", adj = NULL, cex = 0.8)
  text(0.58, 1.045, "TTS+3000", adj = NULL, cex = 0.8)
  text(0.71, 1.045, "TTS+6000", adj = NULL, cex = 0.8)

  # Closing graph and rasterizing image
  dev.off()
  #system(paste("convert -density 300 +antialias -compress lzw -quality 100 a.pdf ",outFileName,sep=""))
  #system("rm a.pdf")

  # Printing clusters to file
  to_write = c()
  groupDict = groups
  names(groupDict) = as.character(groupNames)
  for(i in rev(out$rowInd)){
    rowName = as.character(groupNames[i])
    spl = strsplit(rowName, ":", fixed = TRUE)
    strand = gsub("_NoDuP","",spl[[1]][6])
    to_write = c(to_write, spl[[1]][1], spl[[1]][2], spl[[1]][3], spl[[1]][4], groupDict[rowName], strand)
  }
  to_write = as.character(to_write)
  write(to_write, file = outClusterFileName, ncolumns = 6, append = FALSE, sep = "\t")

}

makeUnique <- function(vector){
  newVector = vector
  while(sum(duplicated(newVector)) != 0){
    dupIndexes = which(duplicated(newVector))
    for(i in 1:length(dupIndexes)){
      newVector[dupIndexes[i]] = paste(newVector[dupIndexes[i]],"NoDuP",sep="_")
    }
  }
  return(newVector)
}

createRowNames <- function(data){
  rowNamesVec = c()
  for(i in 1:nrow(data)){
    chrom = as.character(data[i,1])
    start = as.character(data[i,2])
    end = as.character(data[i,3])
    gene = as.character(data[i,4])
    score = as.character(data[i,5])
    strand = as.character(data[i,6])
    rowNamesVec = c(rowNamesVec, paste(chrom,start,end,gene,score,strand,sep=":"))
  }
  retVec = makeUnique(as.character(rowNamesVec))
  return(retVec)
}

###################################################################################################
# Execution
###################################################################################################

# Reading table
inputTable = read.table(inputTableFileName, sep="\t", header=FALSE)
signalTable = inputTable[,7:ncol(inputTable)]
numberOfElements = nrow(signalTable)

# Make aggregate plot
listOfVectors = vector("list", 1)
signalVector = as.numeric(colMeans(signalTable))
signalVector[is.na(signalVector)] = 0.0
listOfVectors[[1]] = signalVector
yMax = max(signalVector)
yMin = min(signalVector)
linePlot(yMax, yMin, label, listOfVectors, numberOfElements, outputAggrFileName)

# Make heatmap
inputTableNoNA = inputTable[complete.cases(inputTable), ]
rownames(inputTableNoNA) = createRowNames(inputTableNoNA)
signalTable = inputTableNoNA[,7:ncol(inputTableNoNA)]
signalMatrix = as.matrix(signalTable)
rownames(signalMatrix) = rownames(signalTable)
printRowNames = FALSE
if(nrow(signalTable) < 200){
  printRowNames = TRUE
}
createHeatmap(signalMatrix, label, 12, outputHeatFileName, outputClusFileName, printRowNames)
  

