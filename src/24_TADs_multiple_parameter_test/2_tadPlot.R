
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
name = args[3]
treat = args[4]
matrixFileName = args[5]
tadFileName = args[6]
outputFileName = args[7]

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap <- function(matrix, tads, resolution, minValue, name, cond, outputFileName){

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
  par(mar=c(10,10,10,10))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 50,length=101))

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
  text(x = 1.16, y = 1.20, labels = name, col = "black", cex = 2.0)
  text(x = 1.16, y = 1.15, labels = cond, col = "black", cex = 2.0)

  # Black border around plot
  x1 = -0.332
  y1 = 1.332
  x2 = 1.059
  y2 = -0.332
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = "black", lty = 1, lwd = 2)

  # TADs
  for(i in 1:nrow(tads)){
    color = "black"
    bin1 = ((tads[i,1]-minValue-(resolution/2))/resolution)
    bin2 = ((tads[i,2]-minValue-(resolution/2))/resolution)
    createTAD(bin1, bin2, color, 1)
  }

  # Closing colorkey graph
  dev.off()

}

# Creating TAD
createTAD <- function(bin1, bin2, color, mylty){
  #x1 = -0.33167 + (bin1 * 0.00347085)
  #y1 = 1.33167 - (bin1 * 0.0041625)
  #x2 = -0.33167 + (bin2 * 0.00347085)
  #y2 = 1.33167 - (bin2 * 0.0041625)
  x1 = -0.33167 + (bin1 * 0.00346215)
  y1 = 1.33167 - (bin1 * 0.00415205)
  x2 = -0.33167 + (bin2 * 0.00346215)
  y2 = 1.33167 - (bin2 * 0.00415205)
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = color, lty = mylty, lwd = 2)
}

##########################################
### Execution
##########################################

# Reading matrices
matrix = as.matrix(read.table(matrixFileName, sep="\t", header=FALSE))

# Reading tads
tads = read.table(tadFileName, sep="\t", header=TRUE)

# Plotting heatmap
createHeatmap(matrix, tads, resolution, minValue, name, treat, outputFileName)


