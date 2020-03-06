
# Import
import os
import sys

# Table List
il = "/home/egg/Projects/Papantonis_Stag/Results/9_tf_enrichment_new/tables/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/9_tf_enrichment_new/graphs/"

"""
tableList = ["DNase_peaks_filter_1_vs_2", "DNase_summits_filter_1_vs_2", "footprints_1_vs_2", "H3K27ac_peaks_filter_1_vs_2", "H3K27ac_summits_filter_1_vs_2"]
pvalueThresholdList = ["0.01", "0.01", "0.01", "0.01", "0.01"]
graphWidthList = ["6", "6", "6", "6", "6"]
graphHeightList = ["8", "8", "8", "8", "8"]
heatmapTitleList = ["\"TF_EnrichmentNWLDNA_Peaks\"", "\"TF_EnrichmentNWLDNA_Summits\"", "\"TF_EnrichmentNWLDNA_Footprints\"", "\"TF_EnrichmentNWLH3K27ac_Peaks\"", "\"TF_EnrichmentNWLH3K27ac_Summits\""]
heatmapDiffTitleList = ["\"Differential_TF_EnrichmentNWLDNA_Peaks\"", "\"Differential_TF_EnrichmentNWLDNA_Summits\"", "\"Differential_TF_EnrichmentNWLDNA_Footprints\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Peaks\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Summits\""]
label1List = ["\"STAG1\"", "\"STAG1\"", "\"STAG1\"", "\"STAG1\"", "\"STAG1\""]
label2List = ["\"STAG2\"", "\"STAG2\"", "\"STAG2\"", "\"STAG2\"", "\"STAG2\""]
mainTitleSizeList = ["1.0", "1.0", "1.0", "1.0", "1.0"]
keyTitleList = ["\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\""]
keyDiffTitleList = ["\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\""]
keySizeList = ["3", "3", "3", "3", "3"]
rowLabelSizeList = ["0.4", "0.4", "1.1", "0.4", "0.4"]
colLabelSizeList = ["1.0", "1.0", "1.0", "1.0", "1.0"]
xMarginList = ["3", "3", "6.3", "3", "3"]
yMarginList = ["1.6", "1.6", "1.6", "1.6", "1.6"]
lheiXList = ["1.2", "1.2", "1.2", "1.2", "1.2"]
lheiYList = ["6.0", "6.0", "6.0", "6.0", "6.0"]
sepWidthXList = ["0.05", "0.05", "0.05", "0.05", "0.05"]
sepWidthYList = ["0.05", "0.05", "0.05", "0.05", "0.05"]



# Heatmap Lists
tableList = ["DNase_peaks_filter_1_vs_2", "DNase_peaks_filter_random", "DNase_summits_filter_1_vs_2", "DNase_summits_filter_random", "footprints_1_vs_2", "footprints_random", "H3K27ac_peaks_filter_1_vs_2", "H3K27ac_peaks_filter_random", "H3K27ac_summits_filter_1_vs_2", "H3K27ac_summits_filter_random"]
pvalueThresholdList = ["0.001", "0.001", "0.001", "0.001", "0.001", "0.001", "0.001", "0.001", "0.001", "0.001"]
graphWidthList = ["6", "6", "6", "6", "6", "6", "6", "6", "6", "6"]
graphHeightList = ["8", "8", "8", "8", "8", "8", "8", "8", "8", "8"]
heatmapTitleList = ["\"TF_EnrichmentNWLDNA_Peaks\"", "\"TF_EnrichmentNWLDNA_Peaks\"", "\"TF_EnrichmentNWLDNA_Summits\"", "\"TF_EnrichmentNWLDNA_Summits\"", "\"TF_EnrichmentNWLDNA_Footprints\"", "\"TF_EnrichmentNWLDNA_Footprints\"", "\"TF_EnrichmentNWLH3K27ac_Peaks\"", "\"TF_EnrichmentNWLH3K27ac_Peaks\"", "\"TF_EnrichmentNWLH3K27ac_Summits\"", "\"TF_EnrichmentNWLH3K27ac_Summits\""]
heatmapDiffTitleList = ["\"Differential_TF_EnrichmentNWLDNA_Peaks\"", "\"Differential_TF_EnrichmentNWLDNA_Peaks\"", "\"Differential_TF_EnrichmentNWLDNA_Summits\"", "\"Differential_TF_EnrichmentNWLDNA_Summits\"", "\"Differential_TF_EnrichmentNWLDNA_Footprints\"", "\"Differential_TF_EnrichmentNWLDNA_Footprints\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Peaks\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Peaks\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Summits\"", "\"Differential_TF_EnrichmentNWLH3K27ac_Summits\""]
label1List = ["\"STAG1\"", "\"random\"", "\"STAG1\"", "\"random\"", "\"STAG1\"", "\"random\"", "\"STAG1\"", "\"random\"", "\"STAG1\"", "\"random\""]
label2List = ["\"STAG2\"", "\"random\"", "\"STAG2\"", "\"random\"", "\"STAG2\"", "\"random\"", "\"STAG2\"", "\"random\"", "\"STAG2\"", "\"random\""]
mainTitleSizeList = ["0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6", "0.6"]
keyTitleList = ["\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\""]
keyDiffTitleList = ["\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\"", "\"p-value_(-log10)\""]
keySizeList = ["2.5", "2.5", "2.5", "2.5", "2.5", "2.5", "2.5", "2.5", "2.5", "2.5"]
rowLabelSizeList = ["0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4", "0.4"]
colLabelSizeList = ["0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8", "0.8"]
xMarginList = ["3", "3", "3", "3", "3", "3", "3", "3", "3", "3"]
yMarginList = ["3", "3", "3", "3", "3", "3", "3", "3", "3", "3"]
lheiXList = ["1.2", "1.2", "1.2", "1.2", "1.2", "1.2", "1.2", "1.2", "1.2", "1.2"]
lheiYList = ["6.0", "6.0", "6.0", "6.0", "6.0", "6.0", "6.0", "6.0", "6.0", "6.0"]
sepWidthXList = ["0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05"]
sepWidthYList = ["0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05", "0.05"]


tableList = ["DNase_peaks_filter_1_vs_2"]
pvalueThresholdList = ["0.01"]
graphWidthList = ["6"]
graphHeightList = ["8"]
heatmapTitleList = ["\"TF_EnrichmentNWLDNAse_Peaks\"",]
heatmapDiffTitleList = ["\"DifferentialNWLTF_EnrichmentNWLDNAse_Peaks\""]
label1List = ["\"STAG1\""]
label2List = ["\"STAG2\""]
mainTitleSizeList = ["1.0"]
keyTitleList = ["\"p-value_(-log10)\""]
keyDiffTitleList = ["\"p-valueNWL(-log10[STAG1])-(-log10[STAG2])\""]
keySizeList = ["3"]
rowLabelSizeList = ["1.0"]
colLabelSizeList = ["1.0"]
xMarginList = ["4.0"]
yMarginList = ["1.6"]
lheiXList = ["1.0"]
lheiYList = ["6.0"]
sepWidthXList = ["0.05"]
sepWidthYList = ["0.05"]
"""

tableList = ["footprints_1_vs_2"]
pvalueThresholdList = ["0.05"]
graphWidthList = ["6"]
graphHeightList = ["8"]
heatmapTitleList = ["\"TF_EnrichmentNWLDNAse_Footprints\"",]
heatmapDiffTitleList = ["\"DifferentialNWLTF_EnrichmentNWLDNAse_Footprints\""]
label1List = ["\"STAG1\""]
label2List = ["\"STAG2\""]
mainTitleSizeList = ["1.0"]
keyTitleList = ["\"p-value_(-log10)\""]
keyDiffTitleList = ["\"p-valueNWL(-log10[STAG1])-(-log10[STAG2])\""]
keySizeList = ["3"]
rowLabelSizeList = ["0.4"]
colLabelSizeList = ["1.0"]
xMarginList = ["4.0"]
yMarginList = ["1.6"]
lheiXList = ["1.0"]
lheiYList = ["6.0"]
sepWidthXList = ["0.05"]
sepWidthYList = ["0.05"]

# Table Loop
for i in range(0,len(tableList)):

  # Input
  pvalueThreshold = pvalueThresholdList[i]
  graphWidth = graphWidthList[i]
  graphHeight = graphHeightList[i]
  heatmapTitle = heatmapTitleList[i]
  heatmapDiffTitle = heatmapDiffTitleList[i]
  label1 = label1List[i]
  label2 = label2List[i]
  mainTitleSize = mainTitleSizeList[i]
  keyTitle = keyTitleList[i]
  keyDiffTitle = keyDiffTitleList[i]
  keySize = keySizeList[i]
  rowLabelSize = rowLabelSizeList[i]
  colLabelSize = colLabelSizeList[i]
  xMargin = xMarginList[i]
  yMargin = yMarginList[i]
  lheiX = lheiXList[i]
  lheiY = lheiYList[i]
  sepWidthX = sepWidthXList[i]
  sepWidthY = sepWidthYList[i]
  inputFileName = il+tableList[i]+".txt"
  outputFilePrefix = ol+tableList[i]

  # Execution
  command = "R CMD BATCH '--args '"+pvalueThreshold+"' '"+graphWidth+"' '"+graphHeight+"' '"+heatmapTitle+"' '"+heatmapDiffTitle+"' '"+label1+"' '"+label2+"' '"+mainTitleSize+"' '"+keyTitle+"' '"+keyDiffTitle+"' '"+keySize+"' '"+rowLabelSize+"' '"+colLabelSize+"' '"+xMargin+"' '"+yMargin+"' '"+lheiX+"' '"+lheiY+"' '"+sepWidthX+"' '"+sepWidthY+"' '"+inputFileName+"' '"+outputFilePrefix+" 2_heatmap.R 2_heatmap.Rout"
  os.system(command)


