
##########################################
### Initialization
##########################################

# Import
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
args <- commandArgs(trailingOnly = TRUE)
itisstd = args[1]
minV = as.numeric(args[2])
maxV = as.numeric(args[3])
inputTableFileName = args[4]
outputFileNameAll = args[5]
outputFileNameBlue = args[6]
outputFileNameRed = args[7]

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
createHeatmap <- function(minV, maxV, data, outputFile){

  # Initializing graph
  png(filename = outputFile, width = 180 * graphWidth, height = 180 * graphHeight)
  #pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("blue", "white", "red"))(100)
  hmbreaks = c(seq(minV, maxV, length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth, lwid = heatmapLwid, lhei = heatmapLhei, sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo, dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey,
            labRow = heatmapLabRow, labCol = heatmapLabCol)

  colkey(col = colorScheme, clim = c(minV, maxV), clab = NULL, clog = FALSE, add = TRUE, cex.clab = NULL, col.clab = NULL, side.clab = NULL, line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.0, width = 0.8, dist = -0.01, shift = 0, addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE, line = NA, pos = NA, outer = FALSE, font = NA, lty = 2, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL, hadj = 0.0, padj = NA, cex.axis = par("cex.axis"), mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Closing colorkey graph
  dev.off()

}

##########################################
### Execution
##########################################

inputTableAll = as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE))
inputTableBlue = as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE))
inputTableRed = as.matrix(read.table(inputTableFileName, sep="\t", header=FALSE))

if(itisstd == "Y"){
  for(i in 1:nrow(inputTableAll)){
    for(j in 1:ncol(inputTableAll)){
      if(inputTableAll[i,j] == 0){
        inputTableAll[i,j] = 0
        inputTableBlue[i,j] = 0
        inputTableRed[i,j] = 0
      } else if(inputTableAll[i,j] < 0){
        inputTableAll[i,j] = -(10 + log10(-inputTableAll[i,j]))
        inputTableBlue[i,j] = -(10 + log10(-inputTableAll[i,j]))
        inputTableRed[i,j] = 0
      } else if(inputTableAll[i,j] > 0){
        inputTableAll[i,j] = 10 + log10(inputTableAll[i,j])
        inputTableRed[i,j] = 10 + log10(inputTableAll[i,j])
        inputTableBlue[i,j] = 0
      }
    }
  }
}

inputTableBlue[inputTableBlue>0] = 0
inputTableRed[inputTableRed<0] = 0


createHeatmap(minV, maxV, inputTableAll, outputFileNameAll)
createHeatmap(minV, maxV, inputTableBlue, outputFileNameBlue)
createHeatmap(minV, maxV, inputTableRed, outputFileNameRed)


