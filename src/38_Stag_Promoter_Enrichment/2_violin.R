
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
library(ggpubr)
set.seed(111)

###################################################################################################
# Functions
###################################################################################################

# Dodged Barplot
violinPlot <- function(vecX, vecY, yLabel, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY)
  compare_means(vectory ~ vectorx,  data = dataFr)
  my_comparisons <- list( c("CONTROL", "SHARED"), c("CONTROL", "STAG1"), c("CONTROL", "STAG2"), c("SHARED", "STAG1"), c("SHARED", "STAG2"), c("STAG1", "STAG2") )

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, color = vectorx))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + geom_boxplot(width=0.1, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  #pplot = pplot + ylim(0, 1)
  pplot = pplot + xlab("Regions at the promoter of genes")
  pplot = pplot + ylab(yLabel)
  pplot = pplot + stat_compare_means(comparisons = my_comparisons)
  pplot = pplot + theme(legend.position="none", legend.title=element_blank(), axis.text.x = element_text(size=12), 
                        axis.text.y = element_text(size=12), axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), 
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 90, width = 5, height = 5)

}

###################################################################################################
# Execution
###################################################################################################

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/7_Stag_Promoter_Enrichment/1_Promoter_Enrichment_Table/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/7_Stag_Promoter_Enrichment/2_Promoter_Enrichment_Plot/"
nameList1 = c("exp_STAG1_region_STAG1_only", "exp_STAG1_region_STAG1_predominant", "exp_STAG2_region_STAG1_only", "exp_STAG2_region_STAG1_predominant")
nameList2 = c("exp_STAG1_region_STAG2_only", "exp_STAG1_region_STAG2_predominant", "exp_STAG2_region_STAG2_only", "exp_STAG2_region_STAG2_predominant")
nameList3 = c("exp_STAG1_region_shared", "exp_STAG1_region_nonpredominant", "exp_STAG2_region_shared", "exp_STAG2_region_nonpredominant")
nameList4 = c("STAG1_only", "STAG1_predominant", "STAG2_only", "STAG2_predominant")
inputList1 = c()
inputList2 = c()
inputList3 = c()
outputList = c()
for(i in 1:length(nameList1)){
  inputList1 = c(inputList1, paste(il, nameList1[i], ".txt", sep=""))
  inputList2 = c(inputList2, paste(il, nameList2[i], ".txt", sep=""))
  inputList3 = c(inputList3, paste(il, nameList3[i], ".txt", sep=""))
  outputList = c(outputList, paste(ol, nameList4[i], ".pdf", sep=""))
}

for(i in 1:length(inputList1)){

  # Reading tables
  table1 = read.table(inputList1[i], sep="\t", header=TRUE)
  table2 = read.table(inputList2[i], sep="\t", header=TRUE)
  table3 = read.table(inputList3[i], sep="\t", header=TRUE)
  len1 = nrow(table1)
  len2 = nrow(table2)
  len3 = nrow(table3)
  minLen = min(len1, len2, len3)
  if(minLen == len1){
    controlTable = table1
  } else if(minLen == len2){
    controlTable = table2
  } else {
    controlTable = table3
  }

  # X vec
  xVec = c(rep("CONTROL", minLen), rep("SHARED", minLen), rep("STAG1", minLen), rep("STAG2", minLen))

  # Y vec
  contrVec = as.numeric(controlTable[1:minLen,"CONT_EXP"])
  stag1Vec = as.numeric(table1[1:minLen,"STAG_EXP"])
  stag2Vec = as.numeric(table2[1:minLen,"STAG_EXP"])
  sharedVec = as.numeric(table3[1:minLen,"STAG_EXP"])
  yVec = c(contrVec, sharedVec, stag1Vec, stag2Vec)

  # Creating dodged barplot
  outputFileName = outputList[i]
  if(grepl("STAG1_", outputFileName, fixed=TRUE)){
    yLabel = "STAG1 WT/DEG Expression Fold-Change (log2)"
  } else {
    yLabel = "STAG2 WT/DEG Expression Fold-Change (log2)"
  }
  violinPlot(xVec, yVec, yLabel, outputFileName)

}


