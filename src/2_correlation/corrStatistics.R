
rm(list = ls())

# Libraries
library(gplots)
library(RColorBrewer)

# Functions
createBox <- function(matrix, outFileName){
  pdf(file = outFileName, width = 15, height = 6)
  par(mar = c(5,5,4,2), cex.axis = 0.7)
  boxplot(matrix, use.cols = TRUE, xlab=" ", ylab="Correlation per Region")
  grid()
  dev.off()
}

createHeatmap <- function(matrix, outFileName){
  pdf(file = outFileName, width = 8, height = 8)
  par(mar = c(5,5,4,2), cex.axis = 1.0)
  hmbreak = c(seq(0.6,1.0,length=100))
  hmcol = colorRampPalette(brewer.pal(11,"RdBu"))(99)
  heatmap.2(matrix, Rowv = FALSE, Colv = FALSE, dendrogram = 'none', col = hmcol, breaks = hmbreak,
            notecex = 1.0, notecol = "black", trace = "none", cexRow = 1.0, cexCol = 1.0,
            keysize = 1.0, margins = c(4,4), sepwidth = c(2,2), sepcolor = "black", density.info = "none",
            lhei = c(1,5))
  dev.off()
}

# Input
resLoc = "/home/egg/Projects/hic_corr/result_all/"
suf = ".txt"
inList = c("75726", "75728", "75729", "79643", "79644", "79645", "79646", "IMR90")
outLoc = "/home/egg/Projects/hic_corr/graphs/"

# Reading Input Boxplot
resMatrix = matrix(, nrow = 265, ncol = 0)
labels = c()
rowLabels = c()
colLabels = c()
for(i in (1):(length(inList)-1)){
  for(j in (i+1):(length(inList))){
    table = read.table(paste(resLoc,inList[i],"_",inList[j],suf,sep=""), header = FALSE, sep = "\t")
    matrix = table[,2,drop=FALSE]
    resMatrix = cbind(resMatrix,matrix)
    labels = c(labels,paste(inList[i],"\n",inList[j],sep=""))
    rowLabels = table[,1]
    colLabels = c(colLabels,paste(inList[i],"_",inList[j],sep=""))
  }
}
colnames(resMatrix) = labels

# Boxplot
outFileName = paste(outLoc,"boxplot_all.pdf",sep="")
createBox(resMatrix, outFileName)

# Writing Matrix
rownames(resMatrix) = rowLabels
colnames(resMatrix) = colLabels
matrixFileName = paste(outLoc,"all_correlations.txt",sep="")
write.table(resMatrix, file = matrixFileName, sep = "\t", row.names=TRUE, col.names=TRUE)
colnames(resMatrix) = labels

# Reading Input Heatmap
resHeat = matrix(, nrow = length(inList), ncol = length(inList))
colnames(resHeat) = as.character(inList)
rownames(resHeat) = as.character(inList)
for(i in 1:length(inList)){
  for(j in 1:length(inList)){
    if(i == j){
      resHeat[i,j] = 1.0
    }
    else{
      minV = min(inList[i],inList[j])
      maxV = max(inList[i],inList[j])
      resHeat[i,j] = mean(resMatrix[,(paste(minV[1],"\n",maxV[1],sep=""))])
    }
  }
}

# Heatmap
outFileName = paste(outLoc,"heatmap_all.pdf",sep="")
createHeatmap(resHeat, outFileName)

# Writing Matrix
matrixFileName = paste(outLoc,"correlation_heatmap.txt",sep="")
write.table(resHeat, file = matrixFileName, sep = "\t", row.names=TRUE, col.names=TRUE)

rm(list = ls())


