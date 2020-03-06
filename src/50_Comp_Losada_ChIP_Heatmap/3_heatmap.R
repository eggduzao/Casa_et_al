
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
createHeatmap <- function(table, h1, h2, color, outFileName){

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
  colorScheme = colorRampPalette(c("white", color))(100)
  hmbreaks = c(seq(0.0, 1.5,length=101))

  p = levelplot(t(table[nrow(table):1,]), main="Wendt\t\t\t\t\t\t\t\t\t\t\t\tLosada", cex.main=3.0, xlab="", ylab="", scales=list(draw=F), col.regions=colorScheme, cuts=99, at=seq(0,1.5,0.05), aspect="fill",
                panel = function(...) {
                  panel.levelplot(...)
                  panel.abline(v = (ncol(table)/2) + 1, col="black", lwd=3)
                  panel.abline(v = (ncol(table)/4), col="gray", lwd=3, lty=2)
                  panel.abline(v = 3*(ncol(table)/4), col="gray", lwd=3, lty=2)
                  panel.abline(h = h1, col="black", lwd=3)
                  panel.abline(h = h2, col="black", lwd=3)
                }
               )
  print(p)
  dev.off()

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/21_Comp_Losada_ChIP_Heatmap/1_Table/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/21_Comp_Losada_ChIP_Heatmap/3_Plot/"
inFileList = c()
outNameHeatList = c()
nameList = c("W_SA1_L_SA1", "W_SA2_L_SA2") # , "W_SA2_L_SA1", "W_SA2_L_SA2")
colorList = c("dodgerblue4", "darkgreen")
for(name in nameList){
  inFileList = c(inFileList, paste(il, name, "_2.txt", sep=""))
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
  table[table > 1.5] = 1.5
  print(dim(table))

  # Reclaim 200 shift
  table4 = table[,40:250]

  table5 = table[1:nrow(table1),(ncol(table)-10):ncol(table)]
  table6 = table[1:nrow(table1),301:ncol(table)]
  table7 = cbind(table6, table5)

  table8 = table[(nrow(table1)+1):nrow(table),290:ncol(table)]
  table9 = rbind(table7, table8)

  table = cbind(table4, table9)

  # Heatmap Plot
  createHeatmap(table, nrow(table3), nrow(table3) + nrow(table2), colorList[i], outNameHeatList[i])

}



