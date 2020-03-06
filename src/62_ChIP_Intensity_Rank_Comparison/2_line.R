
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecZ, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, fill=vectorz))
  #pplot = pplot + geom_line()
  pplot = pplot + geom_line(data=subset(dataFr,vectorz == "Wendt_Stag1"), color="firebrick3", alpha = 0.75, size = 1)
  pplot = pplot + geom_line(data=subset(dataFr,vectorz == "Wendt_Stag2"), color="darkorange3", alpha = 0.75, size = 1)
  pplot = pplot + geom_line(data=subset(dataFr,vectorz == "Losada_Stag1"), color="dodgerblue3", alpha = 0.75, size = 1)
  pplot = pplot + geom_line(data=subset(dataFr,vectorz == "Losada_Stag2"), color="darkorchid3", alpha = 0.75, size = 1)
  pplot = pplot + theme_classic()
  #pplot = pplot + scale_x_continuous(breaks = c(5, 25, 50, 75, 100))
  #pplot = pplot + xlim(1, 200)
  pplot = pplot + xlab("ChIP-seq peaks")
  pplot = pplot + ylab("ChIP-seq signal intensity")
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=10), 
                        axis.text.y = element_text(size=14), axis.title=element_text(size=16),
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  #pplot = pplot + scale_y_continuous(minor_breaks = seq(0, 100, 5), breaks = seq(0, 100, 10))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 5, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/33_ChIP_Intensity_Rank_Comparison/1_Tables/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/33_ChIP_Intensity_Rank_Comparison/2_Plots/"

# Reading table
tablew1 = read.table(paste(il, "STAG1_full_peaks.txt", sep=""), header = FALSE)
tablew2 = read.table(paste(il, "STAG2_full_peaks.txt", sep=""), header = FALSE)
tablel1 = read.table(paste(il, "SA1_MCF10A_ChIP-seq.txt", sep=""), header = FALSE)
tablel2 = read.table(paste(il, "SA2_MCF10A_ChIP-seq.txt", sep=""), header = FALSE)

# Creating vectors
vectorw1 = sort(log(as.numeric(tablew1[,5])+1), decreasing = TRUE)
vectorw2 = sort(log(as.numeric(tablew2[,5])+1), decreasing = TRUE)
vectorl1 = sort(log(as.numeric(tablel1[,5])+1), decreasing = TRUE)
vectorl2 = sort(log(as.numeric(tablel2[,5])+1), decreasing = TRUE)

# Changing losada
#vectorl1 = vectorl1 + c(2^seq(from = 1.9, to = 0.0, by = -0.019) - 1, rep(0, length(vectorl1)-101))
#vectorl2 = vectorl2 + c(2^seq(from = 1.9, to = 0.0, by = -0.019) - 1, rep(0, length(vectorl2)-101))

# Vector X
vectorX = log(c(1:length(vectorw1), 1:length(vectorw2), 1:length(vectorl1), 1:length(vectorl2)))

# Vector Y
vectorY = c(vectorw1, vectorw2, vectorl1, vectorl2)

# Vector Z
vectorZ = c(rep("Wendt_Stag1",length(vectorw1)), rep("Wendt_Stag2",length(vectorw2)), rep("Losada_Stag1",length(vectorl1)), rep("Losada_Stag2",length(vectorl2)))

# Barplot
outFileName = paste(ol, "plot.pdf", sep="")
barPlot(vectorX, vectorY, vectorZ, outFileName)


