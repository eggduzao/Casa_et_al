
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
args <- commandArgs(trailingOnly = TRUE)
resolution = as.numeric(args[1])
minValue = as.numeric(args[2])
name1 = args[3]
treat1 = args[4]
name2 = args[5]
treat2 = args[6]
matrixFileName1 = args[7]
matrixFileName2 = args[8]
tadFileName1 = args[9]
tadFileName2 = args[10]
diffTadFileName1 = args[11]
diffTadFileName2 = args[12]
outputFileName = args[13]

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap5 <- function(matrix, tads, diffTads, resolution, minValue, maxValue, outputFile){

  # Parameters
  graphWidth = 10
  graphHeight = 10
  heatmapMargins = c(0,10)
  heatmapSepWidth = c(0,0)
  heatmapLwid = c(0.01,10)
  heatmapLhei = c(0.01,10)
  heatmapSepColor = "black"
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "none"
  heatmapRowv = FALSE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapLabRow = FALSE
  heatmapLabCol = FALSE

  # Initializing graph
  postscript(outputFile, width = graphWidth, height = graphHeight, horizontal = FALSE, paper = "special")
  #png(filename = outputFile, width = graphWidth * 100, height = graphHeight * 100)
  #jpeg(filename = outputFile, width = graphWidth * 72, height = graphHeight * 72, quality = 100)
  par(mar=c(10,10,10,10))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 5,length=101))

  # Heatmap
  heatmap.2(matrix, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, dendrogram = heatmapDendrogram,
            Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey, labRow = heatmapLabRow, labCol = heatmapLabCol)

  colkey(col = colorScheme, clim = c(0.0, 5), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL, side = 4, length = 1.4, width = 1.5, dist = 0.06, shift = -0.07,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE, line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 0.2, padj = NA, cex.axis = 1.2,
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  par(xpd=TRUE)

  # Text
  text(x = 1.19, y = 1.28, labels = "STAG2", col = "black", cex = 1.5)
  text(x = 1.19, y = 1.23, labels = "+ AUX", col = "black", cex = 1.5)

  # Black border around plot
  x1 = -0.33035 + (0*0.002777)
  y1 = 1.32995 - (498*0.003330)
  x2 = -0.33035 + (498*0.002777)
  y2 = 1.32995 - (0*0.003330)
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = "black", lty = 1, lwd = 3.0)

  # original TADs
  for(i in 1:nrow(tads)){
    color = "black"
    bin1 = ((tads[i,1]-minValue)/resolution)-1
    bin2 = ((tads[i,2]-minValue)/resolution)-1
    x1 = -0.33035 + (bin1*0.002777)
    y1 = 1.32995 - (bin1*0.003330)
    x2 = -0.33035 + (bin2*0.002777)
    y2 = 1.32995 - (bin2*0.003330)
    rect(x1, y2, x2, y1, angle = 45, col = "NA", border = color, lty = 1, lwd = 3.0)
  }

  # differential TADs
  for(i in 1:nrow(diffTads)){
    color = "green"
    bin1 = ((tads[i,1]-minValue)/resolution)-1
    bin2 = ((tads[i,2]-minValue)/resolution)-1
    x1 = -0.33035 + (bin1*0.002777)
    y1 = 1.32995 - (bin1*0.003330)
    x2 = -0.33035 + (bin2*0.002777)
    y2 = 1.32995 - (bin2*0.003330)
    rect(x1, y2, x2, y1, angle = 45, col = "NA", border = color, lty = 1, lwd = 3.0)
  }

  # Closing colorkey graph
  dev.off()

}

# Regular Heatmap
createHeatmap10 <- function(matrix, tads, diffTads, resolution, minValue, outputFileName){

  # Parameters
  graphWidth = 10
  graphHeight = 10
  heatmapMargins = c(0,10)
  heatmapSepWidth = c(0,0)
  heatmapLwid = c(0.01,10)
  heatmapLhei = c(0.01,10)
  heatmapSepColor = "black"
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "none"
  heatmapRowv = FALSE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapLabRow = FALSE
  heatmapLabCol = FALSE

  # Initializing graph
  pdf(outputFileName, width = graphWidth, height = graphHeight)
  #png(filename = outputFile, width = graphWidth * 90, height = graphHeight * 90)
  #jpeg(filename = outputFile, width = graphWidth * 72, height = graphHeight * 72, quality = 100)
  par(mar=c(10,10,10,10))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 5,length=101))

  # Heatmap
  heatmap.2(matrix, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, dendrogram = heatmapDendrogram,
            Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey, labRow = heatmapLabRow, labCol = heatmapLabCol)

  colkey(col = colorScheme, clim = c(0.0, 5), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL, side = 4, length = 1.4, width = 1.5, dist = 0.06, shift = -0.07,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE, line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 0.2, padj = NA, cex.axis = 1.2,
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  par(xpd=TRUE)

  # Text
  text(x = 1.19, y = 1.28, labels = "STAG2", col = "black", cex = 1.5)
  text(x = 1.19, y = 1.23, labels = "+ AUX", col = "black", cex = 1.5)

  # Black border around plot
  x1 = -0.33
  y1 = 1.33
  x2 = 1.06
  y2 = -0.33
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = "black", lty = 1, lwd = 3.0)

  # TADs
  for(i in 1:nrow(tads)){
    color = "black"
    bin1 = ((tads[i,1]-minValue)/resolution)-1
    bin2 = ((tads[i,2]-minValue)/resolution)-1
    createTAD(bin1, bin2, color)
  }

  # Differential TADs
  for(i in 1:nrow(diffTads)){
    color = "turquoise2"
    bin1 = ((diffTads[i,1]-minValue)/resolution)-1
    bin2 = ((diffTads[i,2]-minValue)/resolution)-1
    createTAD(bin1, bin2, color)
  }

  # Closing colorkey graph
  dev.off()

}

createTAD <- function(bin1, bin2, color){
  x1 = -0.33167 + (bin1 * 0.00138834)
  y1 = 1.33167 - (bin1 * 0.001665)
  x2 = -0.33167 + (bin2 * 0.00138834)
  y2 = 1.33167 - (bin2 * 0.001665)
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = color, lty = 1, lwd = 2.0)
}

standardize <- function(data){
  minV = min(data)
  maxV = max(data)
  newData = (data-minV)/(maxV-minV)
  return(newData)
}

##########################################
### Execution
##########################################

# Reading matrices
matrix1 = as.matrix(read.table(matrixFileName1, sep="\t", header=FALSE))
matrix2 = as.matrix(read.table(matrixFileName2, sep="\t", header=FALSE))

# Reading tads
tads1 = read.table(tadFileName1, sep="\t", header=TRUE)
tads2 = read.table(tadFileName2, sep="\t", header=TRUE)
diffTads1 = read.table(diffTadFileName1, sep="\t", header=TRUE)
diffTads2 = read.table(diffTadFileName2, sep="\t", header=TRUE)

# Plotting heatmap
outputFileName1 = paste(outputFileName, "_1.pdf", sep="")
outputFileName2 = paste(outputFileName, "_2.pdf", sep="")
createHeatmap10(matrix1, tads1, diffTads1, resolution, minValue, outputFileName1)
createHeatmap10(matrix2, tads2, diffTads2, resolution, minValue, outputFileName2)


