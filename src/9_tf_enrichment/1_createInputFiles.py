
# Import
import os
import sys

# DNase List
counter = 1
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf_enrichment_new/input/"
olm = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf_enrichment_new/results_matching/"
ole = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/tf_enrichment_new/results_enrichment/"
foregroundList = ["footprints_stag1.bed", "footprints_stag1.bed", "footprints_stag2.bed", "footprints_stag2.bed",   
                  "DNase_peaks_filter_stag1.bed", "DNase_peaks_filter_stag1.bed", "DNase_peaks_filter_stag2.bed", "DNase_peaks_filter_stag2.bed", 
                  "DNase_summits_filter_stag1.bed", "DNase_summits_filter_stag1.bed", "DNase_summits_filter_stag2.bed", "DNase_summits_filter_stag2.bed", 
                  "H3K27ac_peaks_filter_stag1.bed", "H3K27ac_peaks_filter_stag1.bed", "H3K27ac_peaks_filter_stag2.bed", "H3K27ac_peaks_filter_stag2.bed", 
                  "H3K27ac_summits_filter_stag1.bed", "H3K27ac_summits_filter_stag1.bed", "H3K27ac_summits_filter_stag2.bed", "H3K27ac_summits_filter_stag2.bed"]
backgroundList = ["footprints_stag2.bed", ".", "footprints_stag1.bed", ".",   
                  "DNase_peaks_filter_stag2.bed", ".", "DNase_peaks_filter_stag1.bed", ".", 
                  "DNase_summits_filter_stag2.bed", ".", "DNase_summits_filter_stag1.bed", ".", 
                  "H3K27ac_peaks_filter_stag2.bed", ".", "H3K27ac_peaks_filter_stag1.bed", ".", 
                  "H3K27ac_summits_filter_stag2.bed", ".", "H3K27ac_summits_filter_stag1.bed", "."]
randPropList = [".", "10", ".", "10", 
                ".", "10", ".", "10", 
                ".", "10", ".", "10", 
                ".", "10", ".", "10", 
                ".", "10", ".", "10"]
outLocation = "./input/"

# DNase loop
for i in range(0,len(foregroundList)):

  # Name
  nameF = foregroundList[i].split(".")[0]
  if(backgroundList[i] == "."): nameB = "random"
  else: nameB = backgroundList[i].split(".")[0]
  name = "__".join([nameF, nameB])

  # Input
  fprThreshold = "0.0001"
  pseudocounts = "1"
  multipleTestAlpha = "0.5"
  randProportion = randPropList[i]
  inputFileName = il+foregroundList[i]
  if(backgroundList[i] == "."): backgroundFileName = "."
  else: backgroundFileName = il+backgroundList[i]
  outputLocationMatching = olm+name+"/"
  outputLocationEnrichment = ole+name+"/"

  # Execution
  inFile = open(outLocation+str(counter)+".txt","w")
  inFile.write("\n".join([fprThreshold, pseudocounts, multipleTestAlpha, randProportion, inputFileName, backgroundFileName, outputLocationMatching, outputLocationEnrichment]))
  inFile.close()
  counter += 1


