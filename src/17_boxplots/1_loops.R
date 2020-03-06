
###################################################################################################
# Import
###################################################################################################

# Import
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plot3D)

# Input 1
#plotLabel = "Average"
#prefix1T = "LOOP1_AVERAGE_"
#prefix2T = "LOOP2_AVERAGE_"
#prefix1C = "LOOP1_AVERAGE_CONTROL_"
#prefix2C = "LOOP2_AVERAGE_CONTROL_"
#pprefix = "AVERAGE_PEAK_"
#suffix = "_NORMALIZED_COUNT.RPKM."
#ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_intersection_signal_table.csv"
#mfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_minusOnly_signal_table.csv"
#pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_plusOnly_signal_table.csv"
#ifnp = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_intersection_peak_signal_table.csv"
#outLoc = "/home/egg/Projects/Wendt_Stag/Previous_Results/17_Boxplots/1_loops/stag1_average_signal/"

# Input 2
#plotLabel = "Total"
#prefix1T = "LOOP1_TOTAL_"
#prefix2T = "LOOP2_TOTAL_"
#prefix1C = "LOOP1_TOTAL_CONTROL_"
#prefix2C = "LOOP2_TOTAL_CONTROL_"
#pprefix = "TOTAL_PEAK_"
#suffix = "_NORMALIZED_COUNT.RPM."
#ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_intersection_signal_table.csv"
#mfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_minusOnly_signal_table.csv"
#pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_plusOnly_signal_table.csv"
#ifnp = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG1_RM_STAG1_intersection_peak_signal_table.csv"
#outLoc = "/home/egg/Projects/Wendt_Stag/Previous_Results/17_Boxplots/1_loops/stag1_total_signal/"

# Input 3
#plotLabel = "Average"
#prefix1T = "LOOP1_AVERAGE_"
#prefix2T = "LOOP2_AVERAGE_"
#prefix1C = "LOOP1_AVERAGE_CONTROL_"
#prefix2C = "LOOP2_AVERAGE_CONTROL_"
#pprefix = "AVERAGE_PEAK_"
#suffix = "_NORMALIZED_COUNT.RPKM."
#ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_intersection_signal_table.csv"
#mfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_minusOnly_signal_table.csv"
#pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_plusOnly_signal_table.csv"
#ifnp = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_plusOnly_peak_signal_table.csv"
#outLoc = "/home/egg/Projects/Wendt_Stag/Previous_Results/17_Boxplots/1_loops/stag2_average_signal/"

# Input 4
plotLabel = "Total"
prefix1T = "LOOP1_TOTAL_"
prefix2T = "LOOP2_TOTAL_"
prefix1C = "LOOP1_TOTAL_CONTROL_"
prefix2C = "LOOP2_TOTAL_CONTROL_"
pprefix = "TOTAL_PEAK_"
suffix = "_NORMALIZED_COUNT.RPM."
ifn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_intersection_signal_table.csv"
mfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_minusOnly_signal_table.csv"
pfn = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_plusOnly_signal_table.csv"
ifnp = "/home/egg/Projects/Wendt_Stag/Previous_Results/15_Big_Table/1_tables/R_STAG2_RM_STAG2_plusOnly_peak_signal_table.csv"
outLoc = "/home/egg/Projects/Wendt_Stag/Previous_Results/17_Boxplots/1_loops/stag2_total_signal/"

###################################################################################################
# Functions
###################################################################################################

# barPlot
barPlot <- function(vecX, vecY, vecZ, plotLabel, outFileName){

  # Parameters
  dataFr = data.frame(vectorx = vecX, vectory = vecY, vectorz = vecZ)

  # Plotting graph
  pplot = ggplot(data=dataFr, aes(x=vectorx, y=vectory, color=vectorz))
  pplot = pplot + geom_violin(trim=FALSE, position = position_dodge(0.72))
  pplot = pplot + geom_boxplot(width=0.1, fill="white", position = position_dodge(0.72), show.legend=FALSE)
  pplot = pplot + theme_classic()
  #pplot = pplot + ylim(-2, 3)
  pplot = pplot + xlab(" ")
  pplot = pplot + ylab(paste(plotLabel," Signal Intensity (log10)",sep=""))
  pplot = pplot + guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))
  pplot = pplot + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size=12), 
                        axis.text.y = element_text(size=12), axis.title=element_text(size=16), 
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
  heatmapLwid = c(0.01,10)
  heatmapLhei = c(0.01,10)
  heatmapSepColor = "black"
  heatmapKeySize = 0
  heatmapTrace = "none"
  heatmapTraceCol = NA
  heatmapDensityInfo = "none"
  heatmapDendrogram = "none"
  heatmapRowv = FALSE
  heatmapColv = FALSE
  heatmapKey = FALSE
  heatmapLabRow = FALSE
  heatmapLabCol = FALSE

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
  text(x = 1.16, y = 1.28, labels = "p-value\n(-log10)")

  text(x = -0.215, y = 1.2, labels = "PERS\n(LOOP1)")
  text(x = 0.015, y = 0.92, labels = "GAIN\n(LOOP1)")
  text(x = 0.245, y = 0.64, labels = "LOST\n(LOOP1)")
  text(x = 0.477, y = 0.36, labels = "PERS\n(LOOP2)")
  text(x = 0.709, y = 0.08, labels = "GAIN\n(LOOP2)")
  text(x = 0.939, y = -0.2, labels = "LOST\n(LOOP2)")

  # Stars
  for(i in 1:nrow(corrMatStar)){
    for(j in 1:ncol(corrMatStar)){
      currX = -0.215 + ((i-1)*0.23)
      currY = 1.2 - ((j-1)*0.28)
      if(corrMatStar[i, j] == 3){text(x = currX, y = currY, labels = "***", col = "black", cex = 3)}
      else if(corrMatStar[i, j] == 2){text(x = currX, y = currY, labels = "**", col = "black", cex = 3)}
      else if(corrMatStar[i, j] == 1){text(x = currX, y = currY, labels = "*", col = "black", cex = 3)}
    }
  }

  # Closing colorkey graph
  dev.off()

}

###################################################################################################
# Execution
###################################################################################################

# Reading tables
ti = read.table(ifn, header = TRUE)
tm = read.table(mfn, header = TRUE)
tp = read.table(pfn, header = TRUE)
tip = read.table(ifnp, header = TRUE)

# Signal List
signalList = c("DNase.seq_UW", "H2AFZ_BROAD", "H3K4me1_USC", "H3K4me2_BROAD", "H3K4me3_UW", "H3K9ac_USC", "H3K9me2_BROAD", "H3K9me3_USC", "H4K20me1_BROAD", "H3K27ac_USC", "H3K27me3_BROAD", "H3K36me3_USC", "H3K79me2_BROAD", "ATF3_HAIB", "CBX3_HAIB", "CEBPB_HAIB", "CTCF_BROAD", "CTCF_HAIB", "CTCF_UW", "EGR1_HAIB", "ELF1_HAIB", "EZH2_BROAD", "EZH2phosphoT487_BROAD", "FOSL1_HAIB", "MAX_HAIB", "POL2RA_USC", "POLR2AphosphoS5_HAIB", "REST_HAIB", "SIN3A_HAIB", "SP1_HAIB", "SRF_HAIB", "TCF7L2_USC", "TEAD4_HAIB", "USF1_HAIB", "YY1_HAIB", "ZBTB33_HAIB", "ZNF274_USC", "RAD21_HAIB", "CHIP001_RAD21", "CHIP003_SMC1", "CHIP005_CTCF", "CHIP007_CTCF")

# Iterating on columns
for(signal in signalList){

  signal1T = paste(prefix1T, signal, suffix, sep="")
  signal2T = paste(prefix2T, signal, suffix, sep="")
  signal1C = paste(prefix1C, signal, suffix, sep="")
  signal2C = paste(prefix2C, signal, suffix, sep="")
  signalPeak = paste(pprefix, signal, suffix, sep="")

  # Fetching signal - Normal Loop 1
  si1 = as.numeric(ti[,signal1T])
  si1 = log10(si1[is.finite(si1)])
  sm1 = as.numeric(tm[,signal1T])
  sm1 = log10(sm1[is.finite(sm1)])
  sp1 = as.numeric(tp[,signal1T])
  sp1 = log10(sp1[is.finite(sp1)])

  # Fetching signal - Normal Loop 2
  si2 = as.numeric(ti[,signal2T])
  si2 = log10(si2[is.finite(si2)])
  sm2 = as.numeric(tm[,signal2T])
  sm2 = log10(sm2[is.finite(sm2)])
  sp2 = as.numeric(tp[,signal2T])
  sp2 = log10(sp2[is.finite(sp2)])

  # Fetching signal - Control Loop 1
  sic1 = as.numeric(ti[,signal1C])
  sic1 = log10(sic1[is.finite(sic1)])-2
  smc1 = as.numeric(tm[,signal1C])
  smc1 = log10(smc1[is.finite(smc1)])-2
  spc1 = as.numeric(tp[,signal1C])
  spc1 = log10(spc1[is.finite(spc1)])-2

  # Fetching signal - Control Loop 2
  sic2 = as.numeric(ti[,signal2C])
  sic2 = log10(sic2[is.finite(sic2)])-2
  smc2 = as.numeric(tm[,signal2C])
  smc2 = log10(smc2[is.finite(smc2)])-2
  spc2 = as.numeric(tp[,signal2C])
  spc2 = log10(spc2[is.finite(spc2)])-2

  # Fetching signal - Peak
  sip = as.numeric(tip[,signalPeak])
  sip = log10(sip[is.finite(sip)])

  # Vector Y
  vectorY = c(rowMeans(cbind(si1, si2)), rowMeans(cbind(sm1, sm2)), rowMeans(cbind(sp1, sp2)),
              rowMeans(cbind(sic1, sic2)), rowMeans(cbind(smc1, smc2)), rowMeans(cbind(spc1, spc2)),
              sip)

  # Vector X
  vectorX = c(rep("Signal inside Stag Loops",length(si1)+length(sm1)+length(sp1)),
              rep("Igg Signal inside Stag Loops",length(sic2)+length(smc2)+length(spc2)),
              rep("Signal inside TF Peaks",length(sip)))

  # Vector Z
  vectorZ = c(rep("Persistent",length(si1)), rep("Lost",length(sp1)), rep("Gained",length(sm1)), 
              rep("Persistent",length(sic2)), rep("Lost",length(spc2)), rep("Gained",length(smc2)),
              rep("Peak",length(sip)))

  # Output File
  outFileName = paste(outLoc, signal, "_box.pdf", sep="")

  # Barplot
  barPlot(vectorX, vectorY, vectorZ, plotLabel, outFileName)

  # Calculating correlation
  firstList = vector("list", 6)
  firstList[[1]] = si1
  firstList[[2]] = sm1
  firstList[[3]] = sp1
  firstList[[4]] = si2
  firstList[[5]] = sm2
  firstList[[6]] = sp2
  secondList = vector("list", 6)
  secondList[[1]] = si1
  secondList[[2]] = sm1
  secondList[[3]] = sp1
  secondList[[4]] = si2
  secondList[[5]] = sm2
  secondList[[6]] = sp2

  # Correlation matrix
  corrMat = matrix(rep(1,length(firstList)*length(secondList)), nrow=length(firstList), ncol=length(secondList), byrow = TRUE)
  corrMatStar = matrix(rep(0,length(firstList)*length(secondList)), nrow=length(firstList), ncol=length(secondList), byrow = TRUE)
  for(i in 1:length(firstList)){
    for(j in 1:length(secondList)){
      t = wilcox.test(firstList[[i]], secondList[[j]], alternative = c("two.sided"), paired = FALSE, correct = TRUE, conf.level = 0.99)
      corrMat[i, j] = -log10(t$p.value)
      if(corrMat[i, j] >= 4){corrMatStar[i, j] = 3}
      else if(corrMat[i, j] >= 3){corrMatStar[i, j] = 2}
      else if(corrMat[i, j] >= 2){corrMatStar[i, j] = 1}
    }
  }

  # Heatmap
  outFileName = paste(outLoc, signal, "_corr.pdf", sep="")
  createHeatmap(corrMat, corrMatStar, outFileName)

}


