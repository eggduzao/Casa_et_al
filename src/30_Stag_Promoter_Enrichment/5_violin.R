
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
violinPlot <- function(vecX, vecY, pValueL, pValueG, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72))
  pplot = pplot + geom_boxplot(width=0.1, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  pplot = pplot + ylim(0, 1)
  pplot = pplot + xlab(paste("-log10(Stag < Control) ~ ",round(pValueL,4)," / -log10(Control > Stag) ~ ",round(pValueG,4),sep=""))
  pplot = pplot + ylab("Expression (RPKM)")
  pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=14), 
                        axis.text.y = element_text(size=14), axis.title.x=element_text(size=10), axis.title.y=element_text(size=16), 
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  #pplot = pplot + scale_y_continuous(minor_breaks = seq(-10 , 10, 0.1), breaks = seq(-10, 10, 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 5, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/4_expression_table/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/30_Stag_Promoter_Enrichment/5_plot/"
stagList = c("STAG1", "STAG2")
inputList = c()
outputList = c()
for(stag in stagList){
  inputList = c(inputList, paste(il, stag, "only_active.txt",sep=""))
  outputList = c(outputList, paste(ol, stag, ".pdf",sep=""))
}

for(i in 1:length(inputList)){

  # Reading tables
  table = read.table(inputList[i], sep="\t", header=TRUE)

  # X vec
  xVec = as.character(table[,"TYPE"])

  # Y vec
  yVec = as.numeric(table[,"EXP"])

  # Calculating p values
  vec1 = as.numeric(table[table[,"TYPE"] == "T", "EXP"])
  vec2 = as.numeric(table[table[,"TYPE"] == "C", "EXP"])
  t1 = wilcox.test(vec1, vec2, alternative = c("less"), paired = FALSE, correct = TRUE, conf.level = 0.95)
  t2 = wilcox.test(vec1, vec2, alternative = c("greater"), paired = FALSE, correct = TRUE, conf.level = 0.95)
  pValueL = -log10(t1$p.value)
  pValueG = -log10(t2$p.value)

  # Creating dodged barplot
  outputFileName = outputList[i]
  violinPlot(xVec, yVec, pValueL, pValueG, outputFileName)

}


