
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
inputTableFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/2_Correlation/correlation.txt"
outputFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/2_Correlation/correlation.pdf"

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
createHeatmap <- function(data, outputFile){

  # Initializing graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white","blue4"))(100)
  hmbreaks = c(seq(0.85,1.0,length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks,
            margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei,
            sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
            dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey,
            labRow = heatmapLabRow, labCol = heatmapLabCol)

  colkey(col = colorScheme, clim = c(0.85, 1.0), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.6, width = 1.5, dist = 0.03, shift = 0,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.2, padj = NA, cex.axis = par("cex.axis"),
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Closing colorkey graph
  dev.off()

}

##########################################
### Execution
##########################################

inputTable = as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE))
createHeatmap(inputTable, outputFileName)


