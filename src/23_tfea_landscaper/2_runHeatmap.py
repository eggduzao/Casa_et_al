
# Import
import os
import sys

# Table List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/23_tfea_landscaper/1_table_tfea/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/23_tfea_landscaper/1_table_tfea/"

# Parameters
tableList = ["intersection_minus", "intersection_plus", "minus_plus", "STAG1_STAG2"]
pvalueThresholdList = ["0.05"] * len(tableList)
graphWidthList = ["6"] * len(tableList)
graphHeightList = ["8"] * len(tableList)
heatmapTitleList = ["\"TF_EnrichmentNWLDNAse_Footprints\""] * len(tableList)
heatmapDiffTitleList = ["\"DifferentialNWLTF_EnrichmentNWLDNAse_Footprints\""] * len(tableList)
label1List = ["\"STAG1\""] * len(tableList)
label2List = ["\"STAG2\""] * len(tableList)
mainTitleSizeList = ["1.0"] * len(tableList)
keyTitleList = ["\"p-value_(-log10)\""] * len(tableList)
keyDiffTitleList = ["\"p-valueNWL(-log10[STAG1])-(-log10[STAG2])\""] * len(tableList)
keySizeList = ["3"] * len(tableList)
rowLabelSizeList = ["1.0"] * len(tableList)
colLabelSizeList = ["1.0"] * len(tableList)
xMarginList = ["4.2"] * len(tableList)
yMarginList = ["1.6"] * len(tableList)
lheiXList = ["1.0"] * len(tableList)
lheiYList = ["6.0"] * len(tableList)
sepWidthXList = ["0.05"] * len(tableList)
sepWidthYList = ["0.05"] * len(tableList)

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


