
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
matrixFileName = args[3]
tadFileName = args[4]
outputFileName = args[5]

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap <- function(matrix, tads, resolution, minValue, outputFile){

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
  tempFileName = paste(outputFile,sep="")
  pdf(tempFileName, width = graphWidth, height = graphHeight)
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

  # Test
  #x = -0.33167 + (999 * 0.00138834)
  #y = 1.33167 - (999 * 0.001665)
  #rect(x, y, x+0.00138, y-0.00167, angle = 45, col = "NA", border = "green", lty = 1, lwd = 0.000001)

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
    bin1 = ((tads[i,1]-minValue)/resolution)-1
    bin2 = ((tads[i,2]-minValue)/resolution)-1
    createTAD(bin1, bin2)
  }

  # Closing colorkey graph
  dev.off()
  #system(paste("convert -density 300 +antialias -compress lzw -quality 100 ",tempFileName," ",outputFile,sep=""))
  #system(paste("rm ",tempFileName,sep=""))

}

createTAD <- function(bin1, bin2){
  x1 = -0.33167 + (bin1 * 0.00138834)
  y1 = 1.33167 - (bin1 * 0.001665)
  x2 = -0.33167 + (bin2 * 0.00138834)
  y2 = 1.33167 - (bin2 * 0.001665)
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = "turquoise2", lty = 1, lwd = 2.0)
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

# Reading matrix
matrix = as.matrix(read.table(matrixFileName, sep="\t", header=FALSE))

# Reading tads
tads = read.table(tadFileName, sep="\t", header=TRUE)

# Plotting heatmap
createHeatmap(matrix, tads, resolution, minValue, outputFileName)


