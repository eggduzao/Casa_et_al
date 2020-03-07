
##########################################
### Initialization
##########################################

# Import
rm(list=ls())
library(reshape2)
library(ggplot2)
library(grid)
library(magick)
library(plot3D)
pdf(NULL)

# Input
args <- commandArgs(trailingOnly = TRUE)
fc = as.numeric(args[1])
inputMatrixFileName = args[2]
outputFileName = args[3]
dir.create(dirname(outputFileName), showWarnings = FALSE)

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
insulationPlot <- function(fc, matrix, outputFileName){

  # Creating input for ggplot2
  matrixMelt = melt(matrix)

  # Verifying if fold-change
  if(fc){
    myBreaks = c(0, 1, 2)
    myColors = c("blue", "white", "red")
  } else {
    myBreaks = c(0, 2)
    myColors = c("white", "red")
  }

  # Geom tile plot
  pplot = ggplot(data = matrixMelt, aes(x=Var1, y=Var2, fill=value))
  pplot = pplot + geom_tile()
  pplot = pplot + theme_void()
  pplot = pplot + scale_fill_gradientn(breaks = myBreaks, colours = myColors, na.value = "grey80")
  pplot = pplot + theme(legend.title = element_blank(), legend.position = "NA")

  # Writing plot to file
  par(mar=c(0,5,0,5))
  pdf(outputFileName, width = 10, height = 10)
  print(pplot, vp=viewport(angle = 137.47, width=0.75, height=0.688, x=.500, y=.0))

  dev.off()

}

##########################################
### Execution
##########################################

# Reading matrix
matrix = read.table(inputMatrixFileName, sep="\t", header=FALSE)
rownames(matrix) = seq(1, nrow(matrix), 1)
colnames(matrix) = seq(1, nrow(matrix), 1)

newMat = as.matrix(matrix)

# Plotting heatmap
insulationPlot(fc, newMat, outputFileName)


