
# Import
import os
import sys

# Table List
il = "/home/egg/Projects/Papantonis_Stag/Results/9_tf_enrichment_new/results_enrichment/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/9_tf_enrichment_new/tables/"
el = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/expression/"
tableList = [
("DNase_peaks_filter_stag1__DNase_peaks_filter_stag2", "DNase_peaks_filter_stag2__DNase_peaks_filter_stag1"), ("DNase_peaks_filter_stag1__random", "DNase_peaks_filter_stag2__random"), 
("DNase_summits_filter_stag1__DNase_summits_filter_stag2", "DNase_summits_filter_stag2__DNase_summits_filter_stag1"), ("DNase_summits_filter_stag1__random", "DNase_summits_filter_stag2__random"), 
("footprints_stag1__footprints_stag2", "footprints_stag2__footprints_stag1"), ("footprints_stag1__random", "footprints_stag2__random"), 
("H3K27ac_peaks_filter_stag1__H3K27ac_peaks_filter_stag2", "H3K27ac_peaks_filter_stag2__H3K27ac_peaks_filter_stag1"), ("H3K27ac_peaks_filter_stag1__random", "H3K27ac_peaks_filter_stag2__random"), ("H3K27ac_summits_filter_stag1__H3K27ac_summits_filter_stag2", "H3K27ac_summits_filter_stag2__H3K27ac_summits_filter_stag1"), ("H3K27ac_summits_filter_stag1__random", "H3K27ac_summits_filter_stag2__random")
]

# Table Loop
for table in tableList:

  # Parameters
  sp = table[0].split("__")
  name1 = "_".join(sp[0].split("_")[:-1])
  name2 = sp[1]
  if(name2 != "random"): name2 = "1_vs_2"
  name = "_".join([name1, name2])
  
  # Aux
  folder1 = "/"+table[0].split("__")[0]+"/"
  folder2 = "/"+table[1].split("__")[0]+"/"

  # Input
  aliasFileName = "/home/egg/rgtdata/hg19/alias_human.txt"
  expressionListFileName1 = el+"STAG1.txt"
  expressionListFileName2 = el+"STAG2.txt"
  inputMeTableName1 = il+table[0]+folder1+"fulltest_statistics.txt"
  inputMeTableName2 = il+table[1]+folder2+"fulltest_statistics.txt"
  outputFileName = ol+name+".txt"

  # Execution
  command = "python 2_createHeatmapTable.py "+" ".join([aliasFileName, expressionListFileName1, expressionListFileName2, inputMeTableName1, inputMeTableName2, outputFileName])
  os.system(command)


