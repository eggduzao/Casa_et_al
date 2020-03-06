
##########################################
### Initialization
##########################################

# Import
library(ggplot2)

# Input
args <- commandArgs(trailingOnly = TRUE)
binSize = as.numeric(args[1])
chromosome = args[2]
inputTableFileName1 = args[3]
inputTableFileName2 = args[4]
outputFileName = args[5]

##########################################
### Functions
##########################################

# lineplot
linePlot <- function(vecX, vecY, vecGroup, chromosome, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorgroup = vecGroup)
  maxY = max(vecY)
  maxX = max(vecX)
  minX = min(vecX)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, group=vectorgroup))
  pplot = pplot + geom_line(aes(color=vectorgroup), size = 0.6)
  pplot = pplot + scale_color_manual(values=c("#629294", "#8e8c8d"))
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("Interaction Separation in cis (Mbp)")
  pplot = pplot + ylab("Normalized Interaction Frequency (log)")
  pplot = pplot + geom_text(x=(maxX+minX)/2, y=maxY, label=chromosome, col = "black", size = 8)
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
  #pplot = pplot + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

##########################################
### Execution
##########################################

# Read table
inputTable1 = read.table(inputTableFileName1, sep="\t", header=FALSE)
inputTable2 = read.table(inputTableFileName2, sep="\t", header=FALSE)

# Binning and preparing vectors
scoreVector1 = c(inputTable1[,2])
scoreVector2 = c(inputTable2[,2])
scoreVector1Binned = c()
for(i in seq(from = 1, to = length(scoreVector1), by = binSize)){
  scoreVector1Binned = c(scoreVector1Binned, sum(scoreVector1[i:min(i+binSize, length(scoreVector1))]) / binSize)
}
scoreVector2Binned = c()
for(i in seq(from = 1, to = length(scoreVector2), by = binSize)){
  scoreVector2Binned = c(scoreVector2Binned, sum(scoreVector2[i:min(i+binSize, length(scoreVector2))]) / binSize)
}

# Fetching vectors
vecX = c(1:length(scoreVector1Binned), 1:length(scoreVector2Binned))
vecY = c(log10(as.numeric(scoreVector1Binned)), log10(as.numeric(scoreVector2Binned)))
vecGroup = c(rep("Control",length(scoreVector1Binned)), rep("SiSA",length(scoreVector2Binned)))

# Plot
linePlot(vecX, vecY, vecGroup, chromosome, outputFileName)


