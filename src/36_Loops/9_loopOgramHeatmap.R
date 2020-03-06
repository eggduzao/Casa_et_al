
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
args <- commandArgs(trailingOnly = TRUE)
inputTableFileName = args[1]
regionsFileName = args[2]
outputFileName = args[3]

# Parameters
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

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createHeatmap <- function(data, nRegions, outputFile){

  # Initializing graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("blue", "white", "red"))(100)
  hmbreaks = c(seq(0.2, 1.2,length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks,
            margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei,
            sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
            dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey,
            labRow = heatmapLabRow, labCol = heatmapLabCol)

  colkey(col = colorScheme, clim = c(0.2, 1.2), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.4, width = 1.5, dist = 0.03, shift = 0,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.2, padj = NA, cex.axis = par("cex.axis"),
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  par(xpd=TRUE)
  text(x = 1.19, y = 1.28, labels = paste("n = ", nRegions, sep=""))

  # Closing colorkey graph
  dev.off()

}

##########################################
### Execution
##########################################

inputTable = log10(as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE)))
regionsTable = read.table(regionsFileName, sep="\t", header=FALSE)
nregions = nrow(regionsTable)
createHeatmap(inputTable, nregions, outputFileName)


