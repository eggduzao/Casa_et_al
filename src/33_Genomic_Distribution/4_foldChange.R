
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
library(scales)
library(plyr)
set.seed(111)

# Input
args <- commandArgs(trailingOnly = TRUE)
inputStag1FileName = args[1]
inputStag2FileName = args[2]
outputRegionFileName = args[3]
outputSignalFileName = args[4]

###################################################################################################
# Functions
###################################################################################################

# Dodged Barplot
dodgedBarPlot <- function(vecX, vecY, yFCVecP, yFCVecN, colorVec, yLab, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorp = yFCVecP, vectorn = yFCVecN)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=reorder(vectorx, -vectory), y=vectory))
  pplot = pplot + geom_bar(stat="identity", position=position_dodge(), aes(fill = reorder(vectorx, -vectory)))
  #pplot = pplot + scale_y_continuous(limits = c(0, 5.5))
  pplot = pplot + geom_text(aes(label=vectorp), vjust = -0.9, color="black", position = position_dodge(0.9), size=3.5)
  pplot = pplot + geom_text(aes(label=vectorn), vjust = 1.5, color="black", position = position_dodge(0.9), size=3.5)
  pplot = pplot + scale_fill_manual(values = colorVec)
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("")
  pplot = pplot + ylab(yLab)
  pplot = pplot + geom_hline(yintercept=0, color="black", size=1)
  pplot = pplot + theme(legend.position="NA", legend.title=element_blank(),
                        axis.text.x = element_text(size=10, angle = 45, hjust = 1, color = colorVec),
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE,
                        legend.key = element_rect(size = 1, color = NA), legend.key.size = unit(0.8, "lines"),
                        legend.text = element_text(size = 8))
  #pplot = pplot + guides(fill=guide_legend(nrow=2))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 8, height = 6)

}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
ts1 = read.table(inputStag1FileName, sep="\t", header=FALSE)
ts2 = read.table(inputStag2FileName, sep="\t", header=FALSE)

################
# Regions
################

# Initialization
xVec = c()
yVec = c()
yFCVecP = c()
yFCVecN = c()

# Fetching STAG vectors
s1Regs = as.character(ts1[,1])
s2Regs = as.character(ts2[,1])
s1Tab = table(s1Regs)
s2Tab = table(s2Regs)
totalCount1 = length(s1Regs)
totalCount2 = length(s2Regs)
for(i in 1:length(s1Tab)){
  label = names(s1Tab)[i]
  value1 = as.numeric(s1Tab[i])
  value2 = as.numeric(s2Tab[i])
  fc = round( log2(((value1/totalCount1)+1) / ((value2/totalCount2)+1)) * 100 , 2)
  xVec = c(xVec, gsub("_", "\n", label))
  yVec = c(yVec, fc)
  if(fc > 0){
    yFCVecP = c(yFCVecP, paste("(",fc,")",sep=""))
    yFCVecN = c(yFCVecN, "")
  }
  if(fc < 0){
    yFCVecN = c(yFCVecN, paste("(",fc,")",sep=""))
    yFCVecP = c(yFCVecP, "")
  }
}

# Creating stacked barplot
colorVec = c("darkgoldenrod4", "palegreen4", "palegreen4", "dodgerblue3", "dodgerblue3", "firebrick", "dimgray", "palegreen4", "firebrick", "palegreen4")
yLab = "Stag Region Fold Change (log2(Stag1/Stag2))"
dodgedBarPlot(xVec, yVec, yFCVecP, yFCVecN, colorVec, yLab, outputRegionFileName)

################
# Signal
################

# Initialization
xVec = c()
yVec = c()
yFCVecP = c()
yFCVecN = c()

# Fetching STAG vectors
s1Regs = as.character(ts1[,1])
s2Regs = as.character(ts2[,1])
s1Tab = table(s1Regs)
s2Tab = table(s2Regs)
totalCount1 = sum(as.numeric(ts1[,2]))
totalCount2 = sum(as.numeric(ts2[,2]))
for(i in 1:length(s1Tab)){
  label = names(s1Tab)[i]
  value1 = sum(as.numeric(ts1[ts1[,1]==label,2]))
  value2 = sum(as.numeric(ts2[ts2[,1]==label,2]))
  fc = round( log2(((value1/totalCount1)+1) / ((value2/totalCount2)+1)) * 100 , 2)
  xVec = c(xVec, gsub("_", "\n", label))
  yVec = c(yVec, fc)
  if(fc > 0){
    yFCVecP = c(yFCVecP, paste("(",fc,")",sep=""))
    yFCVecN = c(yFCVecN, "")
  }
  if(fc < 0){
    yFCVecN = c(yFCVecN, paste("(",fc,")",sep=""))
    yFCVecP = c(yFCVecP, "")
  }
}

# Creating stacked barplot
colorVec = c("darkgoldenrod4", "palegreen4", "palegreen4", "dodgerblue3", "dodgerblue3", "firebrick", "dimgray", "firebrick", "palegreen4", "palegreen4")
yLab = "Stag Signal Fold Change (log2(Stag1/Stag2))"
dodgedBarPlot(xVec, yVec, yFCVecP, yFCVecN, colorVec, yLab, outputSignalFileName)


