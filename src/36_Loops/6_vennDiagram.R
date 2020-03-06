
##########################################
### Initialization
##########################################

# Import
library(VennDiagram)

# Parameters
graphWidth = 800
graphHeight = 800
fill1 = "blue"
fill2 = "yellow"
fill3 = "red"
graphAlpha = 0.3
graphScaled = TRUE
printMode = c("raw", "percent")
marginLen = 0.1
sigDigits = 2
graphRes = 90
imageType = "png"

##########################################
### Heatmap Functions
##########################################

# Regular Heatmap
vennDiagram <- function(xVec, yVec, zVec, outputFile){

  # Venn diagram
  venn.diagram(x = list("STAG1\nOnly" = xVec[1]:xVec[2], "STAG2\nOnly" = yVec[1]:yVec[2], "Loops do not overlap Stags" = zVec[1]:zVec[2]), lwd = "0.01",
               height = graphHeight, width = graphWidth, fill = c(fill1, fill2, fill3), cex = 1.5, cat.cex = 1.5, cat.dist = c(0.08, 0.08, 0.03),
               alpha = graphAlpha, scaled = graphScaled, print.mode = printMode, margin = marginLen, cat.pos = c(300, 60, 180),
               sigdigs = sigDigits, filename = outputFile, resolution = graphRes, imagetype = imageType)

  # Remove log file
  system(paste("rm ",outputFile,"*.log",sep=""))

}

##########################################
### Execution
##########################################

# Input
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/5_Loops/5_ChoosingBestLoops/Venn/"
regionList = c("STAG1_minus_auxin_50Kb_5Mb", "STAG1_minus_auxin_50Kb_10Mb", "STAG1_minus_auxin_150Kb_5Mb" , "STAG1_minus_auxin_150Kb_10Mb",
               "STAG1_plus_auxin_50Kb_5Mb", "STAG1_plus_auxin_50Kb_10Mb", "STAG1_plus_auxin_150Kb_5Mb" , "STAG1_plus_auxin_150Kb_10Mb",
               "STAG2_minus_auxin_50Kb_5Mb", "STAG2_minus_auxin_50Kb_10Mb", "STAG2_minus_auxin_150Kb_5Mb" , "STAG2_minus_auxin_150Kb_10Mb",
               "STAG2_plus_auxin_50Kb_5Mb", "STAG2_plus_auxin_50Kb_10Mb", "STAG2_plus_auxin_150Kb_5Mb" , "STAG2_plus_auxin_150Kb_10Mb")
inList = c()
outList = c()
for(region in regionList){
  inList = c(inList, paste(il,region,".txt",sep=""))
  outList = c(outList, paste(il,region,".png",sep=""))
}

for(i in 1:length(inList)){
  inputTableFileName = inList[i]
  outputFileName = outList[i]
  inputTable = read.table(inputTableFileName, sep="\t", header=TRUE)
  xVec = c(inputTable[1,"STAG1"],inputTable[2,"STAG1"])
  yVec = c(inputTable[1,"STAG2"],inputTable[2,"STAG2"])
  zVec = c(inputTable[1,"NONE"],inputTable[2,"NONE"])
  vennDiagram(xVec, yVec, zVec, outputFileName)
}


