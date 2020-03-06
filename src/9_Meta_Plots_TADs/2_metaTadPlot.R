
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
#args <- commandArgs(trailingOnly = TRUE)
#inputTableFileName = args[1]
#outputFileName = args[2]

inputTableFileName = "/home/egg/Desktop/metatad.txt"
outputFileName = "/home/egg/Desktop/metatad.pdf"

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap <- function(data, outputFile){

  # Heatmap Parameters
  graphWidth = 5
  graphHeight = 5
  heatmapMargins = c(0,5)
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
  heatmapNaCol = "gray"

  # Minimum and maximum values
  minV = as.numeric(quantile(data, c(.05), na.rm = TRUE))
  maxV = as.numeric(quantile(data, c(.95), na.rm = TRUE))

  # Initializing graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white","red"))(100)
  hmbreaks = c(seq(minV,maxV,length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth, lwid = heatmapLwid, lhei = heatmapLhei,
            sepcolor = heatmapSepColor, keysize = heatmapKeySize, trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
            dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey, labRow = heatmapLabRow, labCol = heatmapLabCol, 
            na.color=heatmapNaCol)

  colkey(col = colorScheme, clim = c(minV, maxV), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL, side = 4, length = 1.6, width = 1.5, dist = 0.03, shift = 0, line = NA, pos = NA,
         outer = FALSE, font = NA, lty = 1, lwd = 1, lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 0.2, padj = NA,
         cex.axis = par("cex.axis"), mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Closing colorkey graph
  dev.off()

}

##########################################
### Execution
##########################################

inputTable = log10(as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE)))
createHeatmap(inputTable, outputFileName)


