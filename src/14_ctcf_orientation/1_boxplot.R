
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)

# Input Broad
ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/R_STAG1_RM_STAG1_intersection__CTCF_BROAD.txt"
pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127plus.txt:STAG1_regions__CTCF_BROAD.txt"
nfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127-.txt:STAG1_regions__CTCF_BROAD.txt"
outFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_BROAD.pdf"
outTableFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_BROAD.csv"

# Input HAIB
#ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/R_STAG1_RM_STAG1_intersection__CTCF_HAIB.txt"
#pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127plus.txt:STAG1_regions__CTCF_HAIB.txt"
#nfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127-.txt:STAG1_regions__CTCF_HAIB.txt"
#outFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_HAIB.pdf"
#outTableFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_HAIB.csv"

# Input UW
#ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/R_STAG1_RM_STAG1_intersection__CTCF_UW.txt"
#pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127plus.txt:STAG1_regions__CTCF_UW.txt"
#nfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/69_127-.txt:STAG1_regions__CTCF_UW.txt"
#outFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_UW.pdf"
#outTableFileName = "/home/egg/Projects/Wendt_Stag/Previous_Results/14_Ctcf_Orientation/STAG1_CTCF_UW.csv"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecFill, vectorYLabel, vecRawNumbers, outFileName, outTableFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorfill = vecFill)
  dataFr$vectorfill <- factor(dataFr$vectorfill, levels = c("convergent", "divergent", "+ in tandem", "- in tandem", "+ downstream", "+ upstream", "- downstream", "- upstream", "no CTCF"))
  write.csv(dataFr, file = outTableFileName)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, fill=vectorfill))
  pplot = pplot + geom_bar(stat="identity", position=position_dodge())
  pplot = pplot + scale_x_discrete(limits = c("Persistent", "Plus Only", "Minus Only"))
  #pplot = pplot + scale_y_continuous(limits = c(-1, 30))
  pplot = pplot + geom_text(aes(label=vectorYLabel), vjust = -2.0, color="black", position = position_dodge(0.9), size=2.5)
  pplot = pplot + geom_text(aes(label=vecRawNumbers), vjust = -0.5, color="grey30", position = position_dodge(0.9), size=2.5)
  #pplot = pplot + geom_text(x=1.5, y=29, label = "Percentage of Intersections", color="black", size = 3.5)
  #pplot = pplot + geom_text(x=1.5, y=28, label = "Number of Intersections", color="dodgerblue3", size = 3.5)
  pplot = pplot + scale_fill_brewer(palette = "Paired")
  pplot = pplot + theme_classic()
  pplot = pplot + xlab("STAG Occurence")
  pplot = pplot + ylab("Frequency (%)")
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title=element_text(size=16))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 12, height = 6)

}

# Calculates the number of CTCF positions in a table
#NA NA no CTCF
#NA n - upstream
#n NA - downstream
#NA p + upstream
#p NA + downstream
#n p convergent
#p n divergent
#n n - parallel
#p p + parallel
calculateCrcfOrientation <- function(table){
  retVec = c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  dVec = as.character(table[,8])
  uVec = as.character(table[,15])
  for(i in 1:length(dVec)){
    if(is.na(dVec[i]) & is.na(uVec[i])){
      retVec[1] = retVec[1] + 1
    }
    else if(is.na(dVec[i]) & uVec[i] == "N"){
      retVec[2] = retVec[2] + 1
    }
    else if(dVec[i] == "N" & is.na(uVec[i])){
      retVec[3] = retVec[3] + 1
    }
    else if(is.na(dVec[i]) & uVec[i] == "P"){
      retVec[4] = retVec[4] + 1
    }
    else if(dVec[i] == "P" & is.na(uVec[i])){
      retVec[5] = retVec[5] + 1
    }
    else if(dVec[i] == "N" & uVec[i] == "P"){
      retVec[6] = retVec[6] + 1
    }
    else if(dVec[i] == "P" & uVec[i] == "N"){
      retVec[7] = retVec[7] + 1
    }
    else if(dVec[i] == "N" & uVec[i] == "N"){
      retVec[8] = retVec[8] + 1
    }
    else if(dVec[i] == "P" & uVec[i] == "P"){
      retVec[9] = retVec[9] + 1
    }    
  }
  return(retVec)
}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
ift = read.table(ifn, header = TRUE)
pft = read.table(pfn, header = TRUE)
nft = read.table(nfn, header = TRUE)

# Numerical vectors
iVec = as.numeric(calculateCrcfOrientation(ift))
pVec = as.numeric(calculateCrcfOrientation(pft))
nVec = as.numeric(calculateCrcfOrientation(nft))

# Vector X
vectorX = c(rep("Persistent",9), rep("Plus Only",9), rep("Minus Only",9))

# Vector Y
vectorY = c(round((100*iVec)/sum(iVec),2), round((100*pVec)/sum(pVec),2), round((100*nVec)/sum(nVec),2))

# Vector Fill
vectorFill = c(rep(c("no CTCF", "- downstream", "- upstream", "+ downstream", "+ upstream", "divergent", "convergent", "- in tandem", "+ in tandem"),3))

# Y labels
vectorYLabel = paste(as.character(vectorY),"%",sep="")

# Raw numbers
vecRawNumbers = paste("(",as.character(c(iVec, pVec, nVec)),")",sep="") 

# Bar plot
barPlot(vectorX, vectorY, vectorFill, vectorYLabel, vecRawNumbers, outFileName, outTableFileName)


