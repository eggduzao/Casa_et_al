
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
foldChange = as.numeric(args[1])
regionName = gsub("_", " ", args[2])
signalName = gsub("_", " ", args[3])
inputMatrixFileName = args[4]
outputFileName = args[5]
dir.create(dirname(outputFileName), showWarnings = FALSE)

if(foldChange == 0){
  fc = FALSE
} else {
  fc = TRUE
}

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
insulationPlot <- function(fc, regionName, signalName, matrix, outputFileName){

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
  #pplot = pplot + geom_segment(x = 11, y = 11, xend = 1, yend = 21, linetype="dashed", color = "green", size = 1)
  #pplot = pplot + geom_segment(x = 21, y = 21, xend = 1, yend = 41, linetype="dashed", color = "green", size = 1)
  #pplot = pplot + geom_text(x = -6, y = 14, label = paste("Region: ", regionName, sep=""), size = 6, colour = "black", angle = 45)
  #pplot = pplot + geom_text(x = -6, y = 12, label = paste("Signal: ", signalName, sep=""), size = 6, colour = "black", angle = 45)
  #pplot = pplot + geom_text(x = -11, y = 30, label = paste("Region: ", regionName, sep=""), size = 5, colour = "black", angle = 45)
  #pplot = pplot + geom_text(x = -11, y = 27, label = paste("Signal: ", signalName, sep=""), size = 5, colour = "black", angle = 45)
  pplot = pplot + scale_fill_gradientn(breaks = myBreaks, colours = myColors, na.value = "grey80")
  pplot = pplot + theme(legend.title = element_blank(), legend.position = "NA")

  # Writing plot to file
  #ggsave(outputFileName, plot=pplot, device = "pdf", dpi = 90, width = 10, height = 10)
  par(mar=c(0,5,0,5))
  pdf(outputFileName, width = 10, height = 10)
  #jpeg(outputFileName, width = 2000, height = 1800)
  print(pplot, vp=viewport(angle = -45, width=0.75, height=0.688, x=.500, y=.0))

  dev.off()

}

# Treating matrix
treatingMatrix <- function(matrix){

  newMatrix = matrix

  # Making diagonal = NA / bottom = 0 / all other values = log10
  for(i in 1:nrow(newMatrix)){
    for(j in 1:ncol(newMatrix)){
      if(i == j){
        newMatrix[i,j] = NA
      } else if(j < i){
        newMatrix[i,j] = NA
      #} else if(i == j-1){
      #  if(matrix[i,j] > 0){ newMatrix[i,j] = log10(matrix[i,j] / 3) }
      #  else if(matrix[i,j] < 0){ newMatrix[i,j] = -log10(-matrix[i,j] / 3) }
      } else {
        if(matrix[i,j] > 0){ newMatrix[i,j] = log10(matrix[i,j]) }
        else if(matrix[i,j] < 0){ newMatrix[i,j] = -log10(-matrix[i,j]) }
        else if(matrix[i,j] == 0){ newMatrix[i,j] = 0 }
      }
    }
  }

  # Minimum and maximum positive values
  valueVec = as.numeric(newMatrix[newMatrix>0])
  valueVec = valueVec[!is.na(valueVec)]
  minVp = min(valueVec)
  maxVp = max(valueVec)

  # Minimum and maximum negative values
  valueVec = as.numeric(newMatrix[newMatrix<0])
  valueVec = valueVec[!is.na(valueVec)]
  minVn = min(valueVec)
  maxVn = max(valueVec)

  # Making all values between 0 and 1
  for(i in 1:nrow(newMatrix)){
    for(j in 1:ncol(newMatrix)){
      if(i == j){
        newMatrix[i,j] = NA
      } else if(j < i){
        newMatrix[i,j] = NA
      } else if(newMatrix[i,j] > 0){
        newMatrix[i,j] = 2 * ((newMatrix[i,j] - minVp) / (maxVp - minVp))
      } else if(newMatrix[i,j] < 0){
        newMatrix[i,j] = 2 * ((newMatrix[i,j] - minVn) / (maxVn - minVn))
      }
    }
  }

  return(newMatrix)

}

##########################################
### Execution
##########################################

# Reading matrix
matrix = as.matrix(read.table(inputMatrixFileName, sep="\t", header=FALSE))

# Removing lower part of matrix
newMatrix = treatingMatrix(matrix)

# Write new matrix
write(as.numeric(newMatrix), file = paste(outputFileName, ".tsv", sep=""), ncolumns = ncol(newMatrix), append = FALSE, sep = "\t")

# Plotting heatmap
insulationPlot(fc, regionName, signalName, newMatrix, outputFileName)


