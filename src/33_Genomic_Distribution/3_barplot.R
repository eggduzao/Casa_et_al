
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
legendRows = as.numeric(args[1])
inputStag1FileName = args[2]
inputStag2FileName = args[3]
outputBarFileName = args[4]
outputTotalBarFileName = args[5]
outputViolinFileName = args[6]
outputTotalCtcfBarFileName = args[7]
outputCtcfViolinFileName = args[8]

###################################################################################################
# Functions
###################################################################################################

# Stacked Barplot
stackedBarPlot <- function(vecX, vecY, vecZ, legendRows, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Cumsum dataFr
  dataFrS <- arrange(dataFr, vectorx)
  dataFrS <- arrange(dataFrS, vectorz, decreasing = TRUE)
  dataFrCS <- ddply(dataFrS, "vectorx", transform, label_ypos = cumsum(vectory))

  # Plotting graph
  pplot = ggplot(data=dataFrCS, aes(x=vectorx, y=vectory, fill=vectorz))
  pplot = pplot + geom_bar(stat="identity")
  #pplot = pplot + geom_text(aes(y=label_ypos, label=vecYP), vjust=1.5, color="black", size=3.5)
  #pplot = pplot + geom_text(aes(y=label_ypos, label=vecYC), vjust=3.0, color="black", size=3.5)
  pplot = pplot + scale_fill_brewer(palette = "Paired")
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("")
  pplot = pplot + ylab("Overlap with Genomic Regions (%)")
  pplot = pplot + coord_flip()
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10),
                        legend.key = element_rect(size = 1, color = NA), legend.key.size = unit(0.8, "lines"),
                        legend.text = element_text(size = 8))
  pplot = pplot + guides(fill=guide_legend(nrow=legendRows))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

# Dodged Barplot
dodgedBarPlot <- function(vecX, vecY, vecZ, yCVec, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorfill = vecZ)
  #dataFr$vectorfill <- factor(dataFr$vectorfill, levels = c("LOOP", "CTCF", "SMC3", "CTCF+SMC3"))

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, fill=vectorfill))
  pplot = pplot + geom_bar(stat="identity", position=position_dodge())
  pplot = pplot + scale_y_continuous(limits = c(0, 5.5))
  pplot = pplot + geom_text(aes(label=yCVec), vjust = -1.2, color="black", position = position_dodge(0.9), size=2.5)
  pplot = pplot + scale_fill_brewer(palette = "Paired")
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("")
  pplot = pplot + ylab("Total STAG Signal (RPKM, log10)")
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10),
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE,
                        legend.key = element_rect(size = 1, color = NA), legend.key.size = unit(0.8, "lines"),
                        legend.text = element_text(size = 8))
  pplot = pplot + guides(fill=guide_legend(nrow=2))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 8, height = 6)

}

# Violin Plot
violinPlot <- function(vecX, vecY, vecZ, outFileName){

  # Parameters
  library(dplyr)
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, color=vectorz))
  pplot = pplot + geom_boxplot(width=0.6, fill="white", position = position_dodge(0.7), show.legend=TRUE)
  pplot = pplot + theme_classic()
  pplot = pplot + ylim(0.5, 2.5)
  pplot = pplot + xlab("")
  pplot = pplot + ylab("STAG Signal (RPKM, log10)")
  #pplot = pplot + guides(color = guide_legend(nrow=legendRows, override.aes = list(shape = 15, size = 5)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10),
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE,
                        legend.key = element_rect(size = 0.5, color = NA), legend.key.size = unit(0.6, "lines"),
                        legend.text = element_text(size = 8))
  #pplot = pplot + scale_y_continuous(minor_breaks = seq(-10 , 10, 0.1), breaks = seq(-10, 10, 0.5))
  pplot = pplot + guides(fill=guide_legend(nrow=2))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 8, height = 6)

}

# Add alpha to color
addAlpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
ts1 = read.table(inputStag1FileName, sep="\t", header=FALSE)
ts2 = read.table(inputStag2FileName, sep="\t", header=FALSE)

#############################
# STACKED
#############################

# Initialization
xVec = c()
yVec = c()
yPVec = c()
yCVec = c()
zVec = c()
zList1 = c()
zList2 = c()

# Fetching STAG1 vectors
s1Regs = as.character(ts1[,1])
s1Tab = table(s1Regs)
totalCount = length(s1Regs)
for(i in 1:length(s1Tab)){
  label = names(s1Tab)[i]
  total = as.numeric(s1Tab[i])
  xVec = c(xVec, "STAG1")
  yVec = c(yVec, round(100*total/totalCount,2))
  #yPVec = c(yPVec, paste(round(100*total/totalCount,2),"%",sep=""))
  #yCVec = c(yCVec, paste("(",total,")",sep=""))
  perc = paste(round(100*total/totalCount,2),"%",sep="")
  count = paste("(",total,")",sep="")
  zList1 = c(zList1, paste("STAG1: ",perc," ",count,sep=""))
}
names(zList1) = names(s1Tab)

# Fetching STAG2 vectors
s2Regs = as.character(ts2[,1])
s2Tab = table(s2Regs)
totalCount = length(s2Regs)
for(i in 1:length(s2Tab)){
  label = names(s2Tab)[i]
  total = as.numeric(s2Tab[i])
  xVec = c(xVec, "STAG2")
  yVec = c(yVec, round(100*total/totalCount,2))
  #yPVec = c(yPVec, paste(round(100*total/totalCount,2),"%",sep=""))
  #yCVec = c(yCVec, paste("(",total,")",sep=""))
  perc = paste(round(100*total/totalCount,2),"%",sep="")
  count = paste("(",total,")",sep="")
  zList2 = c(zList2, paste("STAG2: ",perc," ",count,sep=""))
}
names(zList2) = names(s2Tab)

# Z vec
for(i in 1:length(s1Tab)){
  label = names(s1Tab)[i]
  percount1 = zList1[i]
  percount2 = zList2[i]
  zVec = c(zVec, paste("  ", gsub("_", " ", label),"  =  ",percount1,"  /  ",percount2,"     ",sep=""))
}

# Creating stacked barplot
stackedBarPlot(xVec, yVec, zVec, legendRows, outputBarFileName)

#############################
# DODGED
#############################

# Initialization
xVec = c()
yVec = c()
yCVec = c()
zVec = c()

# Fetching STAG1 vectors
s1Regs = as.character(ts1[,1])
s1Tab = table(s1Regs)
for(i in 1:length(s1Tab)){
  label = names(s1Tab)[i]
  total = sum(as.numeric(ts1[ts1[,1] == label,2]))
  xVec = c(xVec, "STAG1")
  yVec = c(yVec, log10(total))
  yCVec = c(yCVec, paste("(",round(log10(total),2),")",sep=""))
  zVec = c(zVec, paste(gsub("_", " ", label),"  ",sep=""))
}

# Fetching STAG2 vectors
s2Regs = as.character(ts2[,1])
s2Tab = table(s2Regs)
for(i in 1:length(s2Tab)){
  label = names(s2Tab)[i]
  total = sum(as.numeric(ts2[ts2[,1] == label,2]))
  xVec = c(xVec, "STAG2")
  yVec = c(yVec, log10(total))
  yCVec = c(yCVec, paste("(",round(log10(total),2),")",sep=""))
  zVec = c(zVec, paste(gsub("_", " ", label),"  ",sep=""))
}

# Creating dodged barplot
dodgedBarPlot(xVec, yVec, zVec, yCVec, outputTotalBarFileName)

#############################
# VIOLIN
#############################

# Initialization
xVec = c()
yVec = c()
zVec = c()

# Iterating on STAG1
for(i in 1:nrow(ts1)){
  xVec = c(xVec, "STAG1")
  yVec = c(yVec, log10(as.numeric(ts1[i,2])))
  zVec = c(zVec, paste(gsub("_", " ", ts1[i,1]),"  ",sep=""))
}

# Iterating on STAG2
for(i in 1:nrow(ts2)){
  xVec = c(xVec, "STAG2")
  yVec = c(yVec, log10(as.numeric(ts2[i,2])))
  zVec = c(zVec, paste(gsub("_", " ", ts2[i,1]),"  ",sep=""))
}

# Creating violin plot
violinPlot(xVec, yVec, zVec, outputViolinFileName)


