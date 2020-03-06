
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
width = 20
height = 10
dpi = 90
lowerThreshold = -1
upperThreshold = 1000000000
ylim1 = 5
ylim2 = 7
boxplotWidth = 0.1
xTickSize = 14
xTickAngle = 45
yTickSize = 14
xTitleSize = 1
yTitleSize = 18
minorGridSize = 0.1
majorGridSize = 0.2
xLabel = " "
yLabel = "Loop length (Mbp)"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(table, width, height, dpi, ylim1, ylim2, boxplotWidth, xTickSize, xTickAngle, yTickSize, xTitleSize, yTitleSize, minorGridSize, majorGridSize, minorBreaks, majorBreaks, xLabel, yLabel, outputPdfFileName){

  n = colnames(table)

  # Plotting graph
  pplot = ggplot(data = melt(table), aes(x=variable, y=value))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72), fill = "lightgray")#, aes(fill=variable))
  #pplot = pplot + scale_fill_manual(values = c(n[1] = "black", n[2] = "black", n[3] = "black", n[4] = "black", n[5] = "blue", n[6] = "blue", n[7] = "blue", n[8] = "blue", n[9] = "red", n[10] = "red", n[11] = "red", n[12] = "red", n[13] = "green", n[14] = "green", n[15] = "green", n[16] = "green", n[17] = "orange", n[18] = "orange", n[19] = "orange", n[20] = "orange", n[21] = "purple", n[22] = "purple", n[23] = "purple", n[24] = "purple"))
  #pplot = pplot + scale_fill_manual(values = c("black", "black", "black", "black", "blue", "blue","blue", "blue", "red", "red", "red", "red", "green", "green", "green", "green", "orange", "orange", "orange", "orange", "purple", "purple", "purple", "purple"))
  pplot = pplot + geom_boxplot(width=boxplotWidth, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  pplot = pplot + coord_cartesian(ylim = c(ylim1, ylim2))
  pplot = pplot + xlab(xLabel)
  pplot = pplot + ylab(yLabel)
  pplot = pplot + theme(legend.position="none", axis.text.x = element_text(size=xTickSize, angle = xTickAngle, vjust=1, hjust=1), axis.text.y = element_text(size=yTickSize), 
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
inputFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/24_Wendt_Loop_Length_Distribution/2_Categorized/table.txt"
outputFileName = "/usr/users/egadegu/Projects/Wendt_Stag/Results/24_Wendt_Loop_Length_Distribution/2_Categorized/loop_length_distribution.pdf"

# Reading tables
table = read.table(inputFileName, header = TRUE)
#table2 = table[,c(1,8,15,22, 2,9,16,23, 3,10,17,24, 4,11,18,25, 5,12,19,26, 6,13,20,27, 7,14,21,28)]
table2 = table[,c(1,8,15,22, 2,9,16,23, 4,11,18,25, 5,12,19,26, 6,13,20,27, 7,14,21,28)]

#colnames(table2) = c("SMC3\n(Stag1-AUX)", "SMC3\n(Stag1+AUX)", "SMC3\n(Stag2-AUX)", "SMC3\n(Stag2+AUX)", "RAD21\n(Stag1-AUX)", "RAD21\n(Stag1+AUX)", "RAD21\n(Stag2-AUX)", "RAD21\n(Stag2+AUX)", "CTCF5\n(Stag1-AUX)", "CTCF5\n(Stag1+AUX)", "CTCF5\n(Stag2-AUX)", "CTCF5\n(Stag2+AUX)", "CTCF7\n(Stag1-AUX)", "CTCF7\n(Stag1+AUX)", "CTCF7\n(Stag2-AUX)", "CTCF7\n(Stag2+AUX)", "Shared\n(Stag1-AUX)", "Shared\n(Stag1+AUX)", "Shared\n(Stag2-AUX)", "Shared\n(Stag2+AUX)", "Stag1Only\n(Stag1-AUX)", "Stag1Only\n(Stag1+AUX)", "Stag1Only\n(Stag2-AUX)", "Stag1Only\n(Stag2+AUX)", "Stag2Only\n(Stag1-AUX)", "Stag2Only\n(Stag1+AUX)", "Stag2Only\n(Stag2-AUX)", "Stag2Only\n(Stag2+AUX)")
colnames(table2) = c("SMC3\n(Stag1-AUX)", "SMC3\n(Stag1+AUX)", "SMC3\n(Stag2-AUX)", "SMC3\n(Stag2+AUX)", "RAD21\n(Stag1-AUX)", "RAD21\n(Stag1+AUX)", "RAD21\n(Stag2-AUX)", "RAD21\n(Stag2+AUX)", "CTCF\n(Stag1-AUX)", "CTCF\n(Stag1+AUX)", "CTCF\n(Stag2-AUX)", "CTCF\n(Stag2+AUX)", "Shared\n(Stag1-AUX)", "Shared\n(Stag1+AUX)", "Shared\n(Stag2-AUX)", "Shared\n(Stag2+AUX)", "Stag1Only\n(Stag1-AUX)", "Stag1Only\n(Stag1+AUX)", "Stag1Only\n(Stag2-AUX)", "Stag1Only\n(Stag2+AUX)", "Stag2Only\n(Stag1-AUX)", "Stag2Only\n(Stag1+AUX)", "Stag2Only\n(Stag2-AUX)", "Stag2Only\n(Stag2+AUX)")
table2 = log10(table2) #/ 1000000

# Excluding values based on a threshold
#table[table>upperThreshold] = NA
#table[table<lowerThreshold] = NA

# Barplot
barPlot(table2, width, height, dpi, ylim1, ylim2, boxplotWidth, xTickSize, xTickAngle, yTickSize, xTitleSize, yTitleSize, minorGridSize, majorGridSize, minorBreaks, majorBreaks, xLabel, yLabel, outputFileName)


