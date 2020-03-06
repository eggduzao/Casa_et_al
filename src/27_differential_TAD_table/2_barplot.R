
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

###################################################################################################
# Functions
###################################################################################################

# Dodged Barplot
dodgedBarPlot <- function(vecX, vecY, vecZ, yCVec, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorfill = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=reorder(vectorx, -vectory), y=vectory, fill=vectorfill))
  pplot = pplot + geom_bar(stat="identity", position=position_dodge())
  pplot = pplot + geom_text(aes(label=yCVec), vjust = -0.5, color="black", position = position_dodge(0.9), size=3.5)
  pplot = pplot + scale_fill_brewer(palette = "Paired")
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("")
  pplot = pplot + ylab("Total TADs (%)")
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10),
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE,
                        legend.key = element_rect(size = 1, color = NA), legend.key.size = unit(0.8, "lines"),
                        legend.text = element_text(size = 8))
  pplot = pplot + guides(fill=guide_legend(nrow=2))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 6, height = 5)

}

# Add alpha to color
addAlpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/27_differential_TAD_table/"
condList = c("3B9_5", "69_127", "WT")
binList = c("0_bins", "1_bins", "2_bins", "3_bins", "4_bins", "5_bins")
inputList = c()
outputList = c()
for(cond in condList){
  for(bin in binList){
    inputList = c(inputList, paste(il,bin,"/",cond,".txt",sep=""))
    outputList = c(outputList, paste(il,bin,"/",cond,".pdf",sep=""))
  }
}

for(i in 1:length(inputList)){

  # Reading tables
  table = read.table(inputList[i], sep="\t", header=TRUE)
  countTable = table(table[,6])

  # X vec
  xVec = names(countTable)

  # Y vec
  yVec = as.numeric(countTable)

  # YC vec
  yCVec = as.character(yVec)

  # Y vec
  summ = sum(yVec)
  yVec = 100*yVec/summ
  print(summ)

  # YC vec
  yCVec = paste(yCVec," (",round(yVec,2),"%)",sep="")

  # Z vec
  zVec = paste(gsub("_", " ", xVec),"         ",sep="")

  # X vec
  xVec = gsub("_", "\n", xVec)

  # Creating dodged barplot
  outputFileName = outputList[i]
  dodgedBarPlot(xVec, yVec, zVec, yCVec, outputFileName)

}


