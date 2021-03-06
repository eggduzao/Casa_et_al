
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
# Input
###################################################################################################

# Input
args <- commandArgs(trailingOnly = TRUE)
xLabel = args[1]
yLabel = args[2]
inputFileName = args[3]
outputFileName = args[4]
dir.create(dirname(outputFileName), showWarnings = FALSE)

###################################################################################################
# Functions
###################################################################################################

# Scatter plot
scatterPlot <- function(xVec, yVec, dotCol, xLabel, yLabel, outputFileName){

  # Parameters
  dataFr = data.frame(X=xVec, Y=yVec, Z=dotCol)

  # Calculating correlation
  corrTestSpearman = cor.test(xVec, yVec, alternative = "two.sided", method = "spearman", conf.level = 0.95) # Correlation
  corrSpearman = corrTestSpearman$estimate
  
  # Calculating percentage and total
  totalAA = sum(xVec > 0 & yVec > 0)
  totalAB = sum(xVec > 0 & yVec <= 0)
  totalBA = sum(xVec <= 0 & yVec > 0)
  totalBB = sum(xVec <= 0 & yVec <= 0)
  percAA = paste(round(100*totalAA/length(xVec),2),"%",sep="")
  percAB = paste(round(100*totalAB/length(xVec),2),"%",sep="")
  percBA = paste(round(100*totalBA/length(xVec),2),"%",sep="")
  percBB = paste(round(100*totalBB/length(xVec),2),"%",sep="")
  totalAA = paste("(",totalAA,")",sep="")
  totalAB = paste("(",totalAB,")",sep="")
  totalBA = paste("(",totalBA,")",sep="")
  totalBB = paste("(",totalBB,")",sep="")

  # Plotting graph
  pplot = ggplot(dataFr, aes(x = X, y = Y, color = Z)) 
  pplot = pplot + geom_point(alpha = 0.25)
  pplot = pplot + xlab(xLabel) 
  pplot = pplot + ylab(yLabel)
  pplot = pplot + theme_classic()

  pplot = pplot + geom_hline(yintercept=0, color = "black")
  pplot = pplot + geom_vline(xintercept=0, color = "black")

  pplot = pplot + coord_cartesian(xlim = c(-17, 17), ylim = c(-17, 17))
  pplot = pplot + scale_y_continuous(breaks=seq(-17, 17, 1))
  pplot = pplot + scale_x_continuous(breaks=seq(-17, 17, 1))

  pplot = pplot + annotate("text", label = paste("n = ", length(xVec), sep=""), x = -14.5, y = 16.5, size = 4, colour = "black")
  pplot = pplot + annotate("text", label = paste("r = ", round(corrSpearman, digits = 2), sep=""), x = -15, y = 15.2, size = 4, colour = "black")

  pplot = pplot + annotate("text", label = "AA", x = 16.5, y = 1.0, size = 5, colour = "#00bfc4")
  pplot = pplot + annotate("text", label = "AB", x = 16.5, y = -1.0, size = 5, colour = "#7cae00")
  pplot = pplot + annotate("text", label = "BA", x = -16.5, y = 1.0, size = 5, colour = "#f8766d")
  pplot = pplot + annotate("text", label = "BB", x = -16.5, y = -1.0, size = 5, colour = "#c77cff")

  pplot = pplot + annotate("text", label = percAA, x = 3, y = 16.5, size = 4, colour = "#00bfc4")
  pplot = pplot + annotate("text", label = totalAA, x = 3, y = 15.2, size = 4, colour = "#00bfc4")
  pplot = pplot + annotate("text", label = percAB, x = 3, y = -16.5, size = 4, colour = "#7cae00")
  pplot = pplot + annotate("text", label = totalAB, x = 3, y = -15.2, size = 4, colour = "#7cae00")
  pplot = pplot + annotate("text", label = percBA, x = -3, y = 16.5, size = 4, colour = "#f8766d")
  pplot = pplot + annotate("text", label = totalBA, x = -3, y = 15.2, size = 4, colour = "#f8766d")
  pplot = pplot + annotate("text", label = percBB, x = -3, y = -16.5, size = 4, colour = "#c77cff")
  pplot = pplot + annotate("text", label = totalBB, x = -3, y = -15.2, size = 4, colour = "#c77cff")

  pplot = pplot + theme(plot.title = element_blank(), legend.position="none", axis.ticks.length=unit(-0.25, "cm"), 
                        axis.text.x = element_text(margin=margin(t = 10, r = 10, b = 10, l = 10)), 
                        axis.text.y = element_text(margin=margin(t = 10, r = 10, b = 10, l = 10)),
                        axis.title.x = element_text(margin=margin(t = -4, r = -4, b = -4, l = -4)),
                        axis.title.y = element_text(margin=margin(t = -4, r = -4, b = -4, l = -4)),)

  ggsave(outputFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

addAlpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

###################################################################################################
# Execution
###################################################################################################

# Reading table
table = read.table(inputFileName, sep="\t", header=TRUE)

# Fetching X vector
xVec = as.numeric(table[,"X"])

# Fetching Y vector
yVec = as.numeric(table[,"Y"])

# Colors
zVec = as.numeric(table[,"Z"])
colors = c("firebrick", "darkgoldenrod", "darkgreen", "dodgerblue3")
names(colors) = c(1, 2, 3, 4)

# Fetching dotCol vector
dotCol = colors[zVec]

scatterPlot(xVec, yVec, dotCol, xLabel, yLabel, outputFileName)


