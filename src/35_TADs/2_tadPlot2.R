
##########################################
### Initialization
##########################################

# Import
library(reshape2)
library(ggplot2)

# Input
args <- commandArgs(trailingOnly = TRUE)
resolution = as.numeric(args[1])
minValue = as.numeric(args[2])
name = args[3]
treat = args[4]
matrixFileName = args[5]
tadFileName = args[6]
outputFileName = args[7]

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
createScatterplot <- function(matrix, tads, resolution, minValue, name, cond, outputFileName){

  matrix[matrix > 20] = 20
  melted_cormat = melt(matrix)

  b1 = c()
  b2 = c()
  for(i in 1:nrow(tads)){
    bin1 = ((tads[i,1]-minValue-(resolution/2))/resolution)+1
    bin2 = ((tads[i,2]-minValue-(resolution/2))/resolution)+1
    b1 = c(b1, bin1)
    b2 = c(b2, bin2)
  }

  pplot = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value))
  pplot = pplot + geom_tile()

  pplot = pplot + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                        axis.text.x = element_blank(), axis.text.y = element_blank(),
                        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                        axis.line.x = element_blank(), axis.line.y = element_blank(),
                        legend.title = element_blank(), legend.key.height = unit(2,"cm"), legend.key.width = unit(0.3,"cm"))

  pplot = pplot + scale_fill_gradientn(limits = c(0,20), colours=c("white", "red"))

  print(length(b1))

  pplot = pplot + geom_rect(mapping=aes(xmin=b1[1], xmax=b2[1], ymin=b1[1], ymax=b2[1]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[2], xmax=b2[2], ymin=b1[2], ymax=b2[2]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[3], xmax=b2[3], ymin=b1[3], ymax=b2[3]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[4], xmax=b2[4], ymin=b1[4], ymax=b2[4]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[5], xmax=b2[5], ymin=b1[5], ymax=b2[5]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[6], xmax=b2[6], ymin=b1[6], ymax=b2[6]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[7], xmax=b2[7], ymin=b1[7], ymax=b2[7]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[8], xmax=b2[8], ymin=b1[8], ymax=b2[8]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[9], xmax=b2[9], ymin=b1[9], ymax=b2[9]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[10], xmax=b2[10], ymin=b1[10], ymax=b2[10]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[11], xmax=b2[11], ymin=b1[11], ymax=b2[11]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[12], xmax=b2[12], ymin=b1[12], ymax=b2[12]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[13], xmax=b2[13], ymin=b1[13], ymax=b2[13]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[14], xmax=b2[14], ymin=b1[14], ymax=b2[14]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[15], xmax=b2[15], ymin=b1[15], ymax=b2[15]), fill = NA, color="black", size = 0.3)
  pplot = pplot + geom_rect(mapping=aes(xmin=b1[16], xmax=b2[16], ymin=b1[16], ymax=b2[16]), fill = NA, color="black", size = 0.3)

  #pplot = pplot + geom_rect(mapping=aes(xmin=0, xmax=10, ymin=0, ymax=10), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=10, xmax=20, ymin=10, ymax=20), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=20, xmax=30, ymin=20, ymax=30), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=30, xmax=40, ymin=30, ymax=40), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=40, xmax=50, ymin=40, ymax=50), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=50, xmax=60, ymin=50, ymax=60), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=60, xmax=70, ymin=60, ymax=70), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=70, xmax=80, ymin=70, ymax=80), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=80, xmax=90, ymin=80, ymax=90), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=90, xmax=100, ymin=90, ymax=100), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=100, xmax=110, ymin=100, ymax=110), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=110, xmax=120, ymin=110, ymax=120), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=120, xmax=130, ymin=120, ymax=130), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=130, xmax=140, ymin=130, ymax=140), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=140, xmax=150, ymin=140, ymax=150), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=150, xmax=160, ymin=150, ymax=160), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=160, xmax=170, ymin=160, ymax=170), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=170, xmax=180, ymin=170, ymax=180), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=180, xmax=190, ymin=180, ymax=190), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=190, xmax=200, ymin=190, ymax=200), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=200, xmax=210, ymin=200, ymax=210), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=210, xmax=220, ymin=210, ymax=220), fill = NA, color="black", size = 0.3)
  #pplot = pplot + geom_rect(mapping=aes(xmin=220, xmax=230, ymin=220, ymax=230), fill = NA, color="black", size = 0.3)

  ggsave(outputFileName, plot=pplot, device = "png", dpi = 90, width = 10, height = 10)

}

# Creating TAD
createTAD <- function(bin1, bin2, color, mylty){
  x1 = -0.33167 + (bin1 * 0.00347085)
  y1 = 1.33167 - (bin1 * 0.0041625)
  x2 = -0.33167 + (bin2 * 0.00347085)
  y2 = 1.33167 - (bin2 * 0.0041625)
  rect(x1, y2, x2, y1, angle = 45, col = "NA", border = color, lty = mylty, lwd = 1.5)
}

##########################################
### Execution
##########################################

# Reading matrices
matrix = as.matrix(read.table(matrixFileName, sep="\t", header=FALSE))

# Reading tads
tads = read.table(tadFileName, sep="\t", header=TRUE)

# Plotting heatmap
createScatterplot(matrix, tads, resolution, minValue, name, treat, outputFileName)


