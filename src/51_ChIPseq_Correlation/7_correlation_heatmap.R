
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

###################################################################################################
# Functions
###################################################################################################

# Regular Heatmap
createHeatmap <- function(data, outputFile){

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

  # Initializing graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 1.0,length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth, lwid = heatmapLwid, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey, labRow = heatmapLabRow, labCol = heatmapLabCol)

  # Color key
  colkey(col = colorScheme, clim = c(0, 1), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.4, width = 1.5, dist = 0.03, shift = 0,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.2, padj = NA, cex.axis = par("cex.axis"),
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Col/Row names
  #par(xpd=TRUE)
  #text(x = 1.16, y = 1.28, labels = "p-value\n(-log10)")

  #text(x = -0.215, y = 1.2, labels = "PERS\n(LOOP1)")
  #text(x = 0.015, y = 0.92, labels = "GAIN\n(LOOP1)")
  #text(x = 0.245, y = 0.64, labels = "LOST\n(LOOP1)")
  #text(x = 0.477, y = 0.36, labels = "PERS\n(LOOP2)")
  #text(x = 0.709, y = 0.08, labels = "GAIN\n(LOOP2)")
  #text(x = 0.939, y = -0.2, labels = "LOST\n(LOOP2)")

  # Closing colorkey graph
  dev.off()

}

###################################################################################################
# Execution
###################################################################################################

# Input
inFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/6_Replicate_Correlation/table_stat.txt"
outFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/22_ChIPseq_Correlation/6_Replicate_Correlation/chipseq_replicate_correlation.pdf"

# Reading table
table = as.matrix(read.table(inFileName, header = TRUE, row.names = 1))

# Heatmap
createHeatmap(table, outFileName)


