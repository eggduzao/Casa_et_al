
###################################################################################################
### Parameters
###################################################################################################

# Import
library(gplots)
library(RColorBrewer)
source("/home/egg/R/heatmap.3.R")

# Input
args <- commandArgs(trailingOnly = TRUE)
pvalueThreshold = as.numeric(args[1])
graphWidth = as.numeric(args[2])
graphHeight = as.numeric(args[3])
heatmapTitle = gsub("NWL", "\n", gsub("_", " ", args[4]))
heatmapDiffTitle = gsub("NWL", "\n", gsub("_", " ", args[5]))
label1 = gsub("NWL", "\n", gsub("_", " ", args[6]))
label2 = gsub("NWL", "\n", gsub("_", " ", args[7]))
mainTitleSize = as.numeric(args[8])
keyTitle = gsub("NWL", "\n", gsub("_", " ", args[9]))
keyDiffTitle = gsub("NWL", "\n", gsub("_", " ", args[10]))
keySize = as.numeric(args[11])
rowLabelSize = as.numeric(args[12])
colLabelSize = as.numeric(args[13])
xMargin = as.numeric(args[14])
yMargin = as.numeric(args[15])
lheiX = as.numeric(args[16])
lheiY = as.numeric(args[17])
sepWidthX = as.numeric(args[18])
sepWidthY = as.numeric(args[19])
inputFileName = args[20]
outputFilePrefix = args[21]

###################################################################################################
### Functions
###################################################################################################

# Heatmap dendrogram functions
myDist = function(p1) dist(p1, method="euclidean")
myHclust = function(p2) hclust(p2, method="ward.D")

# Regular Heatmap
createHeatmap <- function(data, colorScheme, hmbreaks, outputFile){

  # Creating graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(cex.main = mainTitleSize)

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks, main = heatmapTitle, margins = c(yMargin, xMargin), dendrogram = "none", offsetRow = -0.2,
            trace = 'none', Rowv = FALSE, Colv = FALSE, density.info = "none", lhei = c(lheiX, lheiY), labCol = c(label1, label2), offsetCol = 0.8,
            sepwidth = c(sepWidthX, sepWidthY), sepcolor = "black", cexCol = colLabelSize, cexRow = rowLabelSize, adjCol = c(0.5,0),
            distfun = myDist, hclustfun = myHclust, keysize = keySize, srtCol = 0, key.title = "NA", key.xlab = keyTitle, key.ylab = "")

  # Closing graph
  dev.off()

}

# Differential Heatmap
createHeatmapDiff <- function(data, colorScheme, hmbreaks, outputFile){

  # Creating graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(cex.main = mainTitleSize)

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks, main = heatmapDiffTitle, margins = c(yMargin-0.5, xMargin), dendrogram = "none",
            trace = 'none', Rowv = FALSE, Colv = FALSE, density.info = "none", lhei = c(lheiX, lheiY), labCol = "", offsetRow = -0.2,
            sepwidth = c(sepWidthX, sepWidthY), sepcolor = "black", cexCol = colLabelSize, cexRow = rowLabelSize, adjCol = c(0.5,0), offsetCol = 0.8,
            distfun = myDist, hclustfun = myHclust, keysize = keySize, srtCol = 0, key.title = "NA", key.xlab = keyDiffTitle, key.ylab = "")

  # Closing graph
  dev.off()

}

###################################################################################################
### Reading input data
###################################################################################################

# Reading Raw Data
table = as.matrix(read.table(inputFileName, sep="\t", header=TRUE, row.names = 1))

# Converting 0's to the minimum value for p-value transformation
#minV = min(table[table>0])
#table[table<0] = minV

###################################################################################################
### Creating heatmap
###################################################################################################

# Keep only the lines in which there is at least one value < p-value1
#filtTable = table[apply(table, 1, function(x) sum(x < pvalueThreshold)) > 0, ]

# Converting all values to log10
#table = -log10(table)
table = cbind(table, table)
table = table[order(table[,1]),]
minV = min(table)
maxV = max(table)

# Parameters
hmcol = colorRampPalette(brewer.pal(11, 'RdBu'))(99)
hmbreaks = c(seq(-0.5, 0.5, length=100))

# Creating heatmaps
outputFileName = paste(outputFilePrefix, ".pdf", sep="")
createHeatmap(table, hmcol, hmbreaks, outputFileName)

###################################################################################################
### Creating differential heatmap
###################################################################################################

# Calulating differences
#diffTable = filtTable[,1] - filtTable[,2]
#diffTable <- cbind(diffTable, diffTable)
#diffTable = diffTable[order(diffTable[,1]),]
#minV = floor(min(diffTable))
#maxV = ceiling(max(diffTable))
#absV = max(abs(minV), abs(maxV))

# Parameters
#hmcol = colorRampPalette(brewer.pal(11, 'RdBu'))(99)
#hmbreaks = c(seq(-absV, absV, length=100))

# Creating Heatmaps
#outputFileName = paste(outputFilePrefix, "_diff.pdf", sep="")
#createHeatmapDiff(diffTable, hmcol, hmbreaks, outputFileName)


