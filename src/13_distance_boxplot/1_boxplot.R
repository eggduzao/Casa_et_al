
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)

# Input K562 DSB
ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/13_Distance_Boxplots/R_STAG1_RM_STAG1_intersection.txt"
mfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/13_Distance_Boxplots/69_127-.txt:STAG1_regions.txt"
pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/13_Distance_Boxplots/69_127plus.txt:STAG1_regions.txt"
outFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/13_Distance_Boxplots/Stag1_Distances.pdf"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory))
  pplot = pplot + geom_violin(trim=FALSE, aes(colour = vectorx))
  pplot = pplot + geom_boxplot(width=0.1, fill="white", show.legend=FALSE)
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("Categories")
  pplot = pplot + ylab("Distance between Loops")
  #pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=14), 
                        axis.text.y = element_text(size=14), axis.title=element_text(size=16))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 6, height = 6)

}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
ifv = as.numeric(read.table(ifn, header = FALSE)[,1])
mfv = as.numeric(read.table(mfn, header = FALSE)[,1])
pfv = as.numeric(read.table(pfn, header = FALSE)[,1])

# Vector X
vectorX = c(rep("Persistent ",length(ifv)), rep("Minus Only ",length(mfv)), rep("Plus Only ",length(pfv)))

# Vector Y
vectorY = c(log10(ifv), log10(mfv), log10(pfv))

barPlot(vectorX, vectorY, outFileName)


