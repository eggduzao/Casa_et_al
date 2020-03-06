
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)
library(reshape2)

# Plot Parameters
width = 5
height = 5
dpi = 90
lowerThreshold = -1
upperThreshold = 1000000000
ylim1 = -20
ylim2 = 100
boxplotWidth = 0.1
xTickSize = 18
xTickAngle = 0
yTickSize = 14
xTitleSize = 18
yTitleSize = 18
minorGridSize = 0.1
majorGridSize = 0.2
xLabel = " "
yLabel = "Average Differential Interactions\n(per chromossome)"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(table, width, height, dpi, ylim1, ylim2, boxplotWidth, xTickSize, xTickAngle, yTickSize, xTitleSize, yTitleSize, minorGridSize, majorGridSize, minorBreaks, majorBreaks, xLabel, yLabel, outputPdfFileName){

  # Plotting graph
  pplot = ggplot(data = melt(table), aes(x=variable, y=value))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72), aes(fill=variable))
  pplot = pplot + geom_boxplot(width=boxplotWidth, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  pplot = pplot + coord_cartesian(ylim = c(ylim1, ylim2))
  pplot = pplot + xlab(xLabel)
  pplot = pplot + ylab(yLabel)
  pplot = pplot + theme(legend.position="none", axis.text.x = element_text(size=xTickSize, angle = xTickAngle), axis.text.y = element_text(size=yTickSize), 
                        axis.title.x=element_text(size=xTitleSize), axis.title.y=element_text(size=yTitleSize), 
                        panel.grid.minor.y = element_line(colour="gray", size=minorGridSize, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=majorGridSize, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  #pplot = pplot + scale_y_continuous(minor_breaks = minorBreaks, breaks = majorBreaks)
  ggsave(outputPdfFileName, plot=pplot, device = "pdf", dpi = dpi, width = width, height = height)

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/usr/users/egadegu/Projects/Wendt_Stag/Results/31_Average_Diff_Interactions_AB/1_Tables/bin/"
ol = "/usr/users/egadegu/Projects/Wendt_Stag/Results/31_Average_Diff_Interactions_AB/2_Plots/bin/"
stagList = c("STAG1+_inter", "STAG1+_intra", "STAG2+_inter", "STAG2+_intra")
ylim1List = c(-125, -270, -5, -110)
ylim2List = c(110, 70, 40, 40)

# Execution
for(i in 1:length(stagList)){

  stag = stagList[i]
  ylim1 = ylim1List[i]
  ylim2 = ylim2List[i]
  inFileName =  paste(il, stag, ".txt", sep="")
  outFileName = paste(ol, stag, ".pdf", sep="")
  
  table = read.table(inFileName, header = TRUE)

  # Barplot
  barPlot(table, width, height, dpi, ylim1, ylim2, boxplotWidth, xTickSize, xTickAngle, yTickSize, xTitleSize, yTitleSize, minorGridSize, majorGridSize, minorBreaks, majorBreaks, xLabel, yLabel, outFileName)

}


