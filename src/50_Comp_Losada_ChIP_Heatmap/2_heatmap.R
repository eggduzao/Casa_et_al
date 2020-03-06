
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

###################################################################################################
# Functions
###################################################################################################

# Heatmap
createHeatmap <- function(table, h1, h2, outFileName){

  # Parameters
  graphWidth = 7
  graphHeight = 10
  heatmapMargins = c(1,1)
  heatmapLmat =  rbind(c(4,3,0), c(2,1,5))
  heatmapLwid = c(0.8, 3.0, 0.2)
  heatmapLhei = c(0.5, 9.5)
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
  heatmapSepColor = c("black", "lightgray", "black")
  heatmapColSep = c(0, 50, 100)
  heatmapSepWidth = c(0.01, 0.01, 0.01)
  heatmapRowSep = c(0, nrow(table))

  # Initializing graph
  png(filename = outFileName, width = graphWidth * 120, height = graphHeight * 120)
  #pdf(outFileName, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 2.0,length=101))

  p = levelplot(t(table[nrow(table):1,]), main="Wendt\t\t\t\t\t\t\t\t\t\t\t\tLosada", cex.main=3.0, xlab="", ylab="", scales=list(draw=F), col.regions=colorScheme, cuts=99, at=seq(0,2,0.05), aspect="fill",
                panel = function(...) {
                  panel.levelplot(...)
                  panel.abline(v = (ncol(table)/2) + 1, col="black", lwd=3)
                  panel.abline(v = (ncol(table)/4), col="dodgerblue", lwd=3, lty=2)
                  panel.abline(v = 3*(ncol(table)/4), col="dodgerblue", lwd=3, lty=2)
                  panel.abline(h = h1, col="black", lwd=3)
                  panel.abline(h = h2, col="black", lwd=3)
                }
)
  print(p)
  # Heatmap
  #heatmap.2(table, col = colorScheme, breaks = hmbreaks, margins = heatmapMargins, sepwidth = heatmapSepWidth, lwid = heatmapLwid, lhei = heatmapLhei,
  #          sepcolor = heatmapSepColor, keysize = heatmapKeySize, trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
  #          dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey, labRow = heatmapLabRow, labCol = heatmapLabCol,
  #          colsep = heatmapColSep, rowsep = heatmapRowSep, lmat = heatmapLmat)

  # Closing graph and rasterizing image
  dev.off()

}

# Scale
specialScale <- function(table){

  minL = 99999
  maxL = -9999
  minR = 99999
  maxR = -9999
  midTable = ncol(table) / 2
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      if(table[i,j] <= midTable){
        if(table[i,j] < minL){minL = table[i,j]}
        if(table[i,j] > maxL){maxL = table[i,j]}
      } else {
        if(table[i,j] < minR){minR = table[i,j]}
        if(table[i,j] > maxR){maxR = table[i,j]}
      }
    }
  }
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      if(table[i,j] <= midTable){
        table[i,j] = (table[i,j] - minL) / (maxL - minL)
      } else {
        table[i,j] = (table[i,j] - minR) / (maxR - minR)
      }
    }
  }

  table = table + min(table)

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/21_Comp_Losada_ChIP_Heatmap/1_Table/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/21_Comp_Losada_ChIP_Heatmap/2_Plot/"
inFileList = c()
outNameHeatList = c()
nameList = c("W_SA1_L_SA1", "W_SA2_L_SA2") # , "W_SA2_L_SA1", "W_SA2_L_SA2")
for(name in nameList){
  inFileList = c(inFileList, paste(il, name, ".txt", sep=""))
  outNameHeatList = c(outNameHeatList, paste(ol, name, ".png", sep=""))
}

# Iterating on input
for(i in 1:length(inFileList)){

  # Reading table
  inputTable = read.table(inFileList[i], header = FALSE, sep = "\t")
  #hist(as.numeric(as.matrix(inputTable[inputTable > 10])))
  table1 = inputTable[inputTable[,1] == 1,2:ncol(inputTable)]
  table2 = inputTable[inputTable[,1] == 2,2:ncol(inputTable)]
  table3 = inputTable[inputTable[,1] == 3,2:ncol(inputTable)]

  # Coercing all values > 10 to 10
  table1[table1 > 10] = 10
  table2[table2 > 10] = 10
  table3[table3 > 10] = 10

  # Sorting table
  table1 = table1[order(rowSums(table1[,c(100:150, 350:400)]),decreasing=T),]
  table2 = table2[order(rowSums(table2[,c(100:150, 350:400)]),decreasing=T),]
  table3 = table3[order(rowSums(table3[,c(100:150, 350:400)]),decreasing=T),]
  table = as.matrix(rbind(table1, table2, table3))
  table[table > 2] = 2

  # Reclaim 200 shift
  table1 = table[,43:250]
  table2 = table[,293:ncol(table)]
  table = as.matrix(cbind(table1, table2))

  # Heatmap Plot
  createHeatmap(table, nrow(table3), nrow(table3) + nrow(table2), outNameHeatList[i])

}



