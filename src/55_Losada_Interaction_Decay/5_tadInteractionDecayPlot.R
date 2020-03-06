
##########################################
### Initialization
##########################################

# Import
library(ggplot2)

# Input
args <- commandArgs(trailingOnly = TRUE)
chromosome = args[1]
resolution = as.numeric(args[2])
inputSA1TableFileName = args[3]
inputSA2TableFileName = args[4]
outputFileName = args[5]

##########################################
### Functions
##########################################

# lineplot
linePlot <- function(vecX, vecY, vecZ, chromosome, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)
  maxY = max(vecY)
  maxX = max(vecX)
  minX = min(vecX)
  mylim1 = -1.2  # -6 -3 -1.2
  mylim2 = 0.5 # 0 0.5 1

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, z=vectorz))
  pplot = pplot + geom_line(aes(color=vectorz), size = 0.6)
  pplot = pplot + scale_color_manual(values=c("#629294", "#8e8c8d"))
  pplot = pplot + theme_classic()
  pplot = pplot + ylim(mylim1, mylim2)
  pplot = pplot + geom_hline(yintercept=0, linetype="dashed", color = "red", size = 0.3, alpha = 0.5)
  #pplot = pplot + xlab("Interaction Separation in cis (Kbp)")
  #pplot = pplot + ylab("Normalized Interaction Frequency (log)")
  pplot = pplot + geom_text(x=(maxX+minX)/2, y=mylim2, label=chromosome, col = "black", size = 8)
  #pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
  pplot = pplot + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

##########################################
### Execution
##########################################

# Read table
sa1Table = read.table(inputSA1TableFileName, sep="\t", header=TRUE)
sa2Table = read.table(inputSA2TableFileName, sep="\t", header=TRUE)

# Vector X
vecX = c(as.numeric(sa1Table[,1]*resolution/1000000), as.numeric(sa2Table[,1]*resolution/1000000))

# Vector Y
vecY = c(log(as.numeric(sa1Table[,3])), log(as.numeric(sa2Table[,3])))

# Vector Z
vecZ = c(rep("STAG1",nrow(sa1Table)), rep("STAG2",nrow(sa2Table)))

# Plot
linePlot(vecX, vecY, vecZ, chromosome, outputFileName)


