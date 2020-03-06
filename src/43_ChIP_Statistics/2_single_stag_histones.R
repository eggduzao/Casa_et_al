
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input
is1 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG1only_mrg_replicates_signal_table.csv"
is2 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG2only_mtg_replicates_signal_table.csv"
isp = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG1only_mrg_replicates_peak_signal_table.csv"
outLoc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/17_boxplots/2_single_stag/average_signal/"
#is1 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG1bestpeaks_filter_signal_table.csv"
#is2 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG2bestpeaks_filter_signal_table.csv"
#isp = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/15_big_statistics_table/1_tables/STAG1bestpeaks_filter_peak_signal_table.csv"
#outLoc = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/17_boxplots/3_single_stag_bestpeaks/total_signal/"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecZ, pValueL, pValueG, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, color=vectorz))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72))
  pplot = pplot + geom_boxplot(width=0.1, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  #pplot = pplot + ylim(-2, 3)
  pplot = pplot + xlab(paste("-log10(SA1 < SA2) ~ ",round(pValueL,4)," / -log10(SA1 > SA2) ~ ",round(pValueG,4),sep=""))
  pplot = pplot + ylab("Average Signal Intensity (log10)")
  pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=14), 
                        axis.text.y = element_text(size=14), axis.title=element_text(size=16), 
                        panel.grid.minor.y = element_line(colour="gray", size=0.2, linetype = "dashed"),
                        panel.grid.major.y = element_line(colour="black", size=0.3, linetype = "dashed"),
                        panel.background = element_rect(fill = NA), panel.ontop = TRUE)
  pplot = pplot + scale_y_continuous(minor_breaks = seq(-10 , 10, 0.1), breaks = seq(-10, 10, 0.5))
  ggsave(outFileName, plot=pplot, device = "pdf", dpi = 300, width = 8, height = 6)

}

# Regular Heatmap
createHeatmap <- function(data, corrMatStar, outputFile){

  # Parameters
  graphWidth = 5
  graphHeight = 5
  heatmapMargins = c(0,5)
  heatmapSepWidth = c(0,0)
  heatmapLwid = c(4,10)
  heatmapLhei = c(1,10)
  heatmapSepColor = "black"
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "none"
  heatmapRowv = FALSE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapLabRow = TRUE
  heatmapLabCol = TRUE

  # Initializing graph
  pdf(outputFile, width = graphWidth, height = graphHeight)
  par(mar=c(5,5,5,5))

  # Color scheme
  colorScheme = colorRampPalette(c("white", "red4"))(100)
  hmbreaks = c(seq(0.0, 5.0,length=101))

  # Heatmap
  heatmap.2(data, col = colorScheme, breaks = hmbreaks,
            margins = heatmapMargins, sepwidth = heatmapSepWidth,
            lwid = heatmapLwid, lhei = heatmapLhei,
            sepcolor = heatmapSepColor, keysize = heatmapKeySize,
            trace = heatmapTrace, tracecol = heatmapTraceCol, density.info = heatmapDensityInfo,
            dendrogram = heatmapDendrogram, Rowv = heatmapRowv, Colv = heatmapColv, key=heatmapKey,
            labRow = heatmapLabRow, labCol = heatmapLabCol)

  # Color key
  colkey(col = colorScheme, clim = c(0, 5), clab = NULL, clog = FALSE, add = TRUE, 
         cex.clab = NULL, col.clab = NULL, side.clab = NULL, 
         line.clab = NULL, adj.clab = NULL, font.clab = NULL,
         side = 4, length = 1.4, width = 1.5, dist = 0.03, shift = 0,
         addlines = FALSE, breaks = NULL, at = NULL, labels = TRUE, tick = TRUE,
         line = NA, pos = NA, outer = FALSE, font = NA, lty = 1, lwd = 1, 
         lwd.ticks = 1, col.axis = NULL, col.ticks = NULL, col.box = NULL,
         hadj = 0.2, padj = NA, cex.axis = par("cex.axis"),
         mgp = NULL, tck = NULL, tcl = NULL, las = 1)

  # Col/Row names
  par(xpd=TRUE)
  text(x = 1.17, y = 1.28, labels = "p-value\n(-log10)")

  for(i in 1:nrow(corrMat)){
    text(x = 0.14, y = 1.158 - ((i-1)*0.041), labels = rownames(data)[i], cex = 0.6, pos = 2)
  }

  text(x = 0.36, y = 1.23, labels = "STAG1 < STAG2", cex = 0.8)
  text(x = 0.82, y = 1.23, labels = "STAG1 > STAG2", cex = 0.8)

  for(i in 1:nrow(corrMat)){
    currY = 1.158 - ((i-1)*0.041)
    if(corrMatStar[i, 1] == 3){text(x = 0.36, y = currY, labels = "***", col = "black", cex = 1)}
    else if(corrMatStar[i, 1] == 2){text(x = 0.36, y = currY, labels = "**", col = "black", cex = 1)}
    else if(corrMatStar[i, 1] == 1){text(x = 0.36, y = currY, labels = "*", col = "black", cex = 1)}
    if(corrMatStar[i, 2] == 3){text(x = 0.82, y = currY, labels = "***", col = "black", cex = 1)}
    else if(corrMatStar[i, 2] == 2){text(x = 0.82, y = currY, labels = "**", col = "black", cex = 1)}
    else if(corrMatStar[i, 2] == 1){text(x = 0.82, y = currY, labels = "*", col = "black", cex = 1)}
  }

  # Closing colorkey graph
  dev.off()

}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
s1 = read.table(is1, header = TRUE)
s2 = read.table(is2, header = TRUE)
sp = read.table(isp, header = TRUE)

# Columns
prefixT = "AVERAGE_"
prefixC = "AVERAGE_CONTROL_"
prefixP = "AVERAGE_PEAK_"
suffix = "_NORMALIZED_COUNT.RPKM."
#prefixT = "TOTAL_"
#prefixC = "TOTAL_CONTROL_"
#prefixP = "TOTAL_PEAK_"
#suffix = "_NORMALIZED_COUNT.RPM."
signalList = c("DNase.seq_UW", "H2AFZ_BROAD", "H3K4me1_USC", "H3K4me2_BROAD", "H3K4me3_UW", "H3K9ac_USC", "H3K9me2_BROAD", "H3K9me3_USC", "H4K20me1_BROAD", "H3K27ac_USC", "H3K27me3_BROAD", "H3K36me3_USC", "H3K79me2_BROAD", "ATF3_HAIB", "CBX3_HAIB", "CEBPB_HAIB", "CTCF_BROAD", "CTCF_HAIB", "CTCF_UW", "EGR1_HAIB", "ELF1_HAIB", "EZH2_BROAD", "EZH2phosphoT487_BROAD", "FOSL1_HAIB", "MAX_HAIB", "POL2RA_USC", "POLR2AphosphoS5_HAIB", "REST_HAIB", "SIN3A_HAIB", "SP1_HAIB", "SRF_HAIB", "TCF7L2_USC", "TEAD4_HAIB", "USF1_HAIB", "YY1_HAIB", "ZBTB33_HAIB", "ZNF274_USC")

# Correlation matrix
corrMat = matrix(rep(1,length(signalList)*2), nrow=length(signalList), ncol=2, byrow = TRUE)
#labels1 = gsub("\.","-",signalList)
#labels2 = gsub("_"," (",labels1)
#labels3 = paste(labels2,")",sep="")
rownames(corrMat) = signalList
colnames(corrMat) = c("STAG1 < STAG2", "STAG1 > STAG2")

# Iterating on columns
counter = 1
for(signal in signalList){

  signalNameT = paste(prefixT, signal, suffix, sep="")
  signalNameC = paste(prefixC, signal, suffix, sep="")
  signalNameP = paste(prefixP, signal, suffix, sep="")

  # Fetching signal - Normal
  si1 = as.numeric(s1[,signalNameT])
  si1 = log10(si1[is.finite(si1)])
  si2 = as.numeric(s2[,signalNameT])
  si2 = log10(si2[is.finite(si2)])

  # Fetching signal - Control
  sic1 = as.numeric(s1[,signalNameC])
  sic1 = log10(sic1[is.finite(sic1)])-2
  sic2 = as.numeric(s2[,signalNameC])
  sic2 = log10(sic2[is.finite(sic2)])-2

  # Fetching signal - Peak
  sip = as.numeric(sp[,signalNameP])
  sip = log10(sip[is.finite(sip)])

  # Vector Y
  vectorY = c(si1, sic1, si2, sic2, sip)

  # Vector X
  vectorX = c(rep("STAG1",length(si1)+length(sic1)), rep("STAG2",length(si2)+length(sic2)), rep("PEAK", length(sip)))

  # Vector Z
  vectorZ = c(rep("Signal inside Stag Loops",length(si1)), rep("Igg Signal inside Stag Loops",length(sic1)),
              rep("Signal inside Stag Loops",length(si2)), rep("Igg Signal inside Stag Loops",length(sic2)),
              rep("Signal inside TF Peaks",length(sip)))

  # Calculating correlation
  t1 = wilcox.test(si1, si2, alternative = c("less"), paired = FALSE, correct = TRUE, conf.level = 0.999)
  t2 = wilcox.test(si1, si2, alternative = c("greater"), paired = FALSE, correct = TRUE, conf.level = 0.999)
  pValueL = -log10(t1$p.value)
  pValueG = -log10(t2$p.value)
  corrMat[counter, 1] = pValueL
  corrMat[counter, 2] = pValueG
  counter = counter + 1

  # Output File
  outFileName = paste(outLoc, signal, "_box.pdf", sep="")

  # Barplot
  barPlot(vectorX, vectorY, vectorZ, pValueL, pValueG, outFileName)

}

# Correlation stars
corrMatStar = matrix(rep(0,nrow(corrMat)*ncol(corrMat)), nrow=nrow(corrMat), ncol=ncol(corrMat), byrow = TRUE)
for(i in 1:nrow(corrMat)){
  for(j in 1:ncol(corrMat)){
    if(corrMat[i, j] >= 4){corrMatStar[i, j] = 3}
    else if(corrMat[i, j] >= 3){corrMatStar[i, j] = 2}
    else if(corrMat[i, j] >= 2){corrMatStar[i, j] = 1}
  }
}
outFileName = paste(outLoc, "correlation.pdf", sep="")
createHeatmap(corrMat, corrMatStar, outFileName)


