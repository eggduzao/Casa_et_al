
# Import
import os
import sys

# Hic files
ml = "/home/egg/Projects/Papantonis_Stag/Data/stag_matrix_files/"
rl = "/home/egg/Projects/Papantonis_Stag/Results/differential_loop_o_gram/intersections/"
matrixFileNameList = [ml+"3B9_5-.txt", ml+"3B9_5-.txt", ml+"3B9_5-.txt", ml+"3B9_5plus.txt", ml+"3B9_5plus.txt", ml+"3B9_5plus.txt", ml+"69_127-.txt", ml+"69_127-.txt", ml+"69_127-.txt", ml+"69_127plus.txt", ml+"69_127plus.txt", ml+"69_127plus.txt"]
regionFileNameList = [rl+"STAG1_3B9_intersection.bed", rl+"STAG1_3B9_minusOnly.bed", rl+"STAG1_3B9_plusOnly.bed", rl+"STAG1_3B9_intersection.bed", rl+"STAG1_3B9_minusOnly.bed", rl+"STAG1_3B9_plusOnly.bed", rl+"STAG1_127_intersection.bed", rl+"STAG1_127_minusOnly.bed", rl+"STAG1_127_plusOnly.bed", rl+"STAG1_127_intersection.bed", rl+"STAG1_127_minusOnly.bed", rl+"STAG1_127_plusOnly.bed"]
outputLocM = "/home/egg/Projects/Papantonis_Stag/Results/differential_loop_o_gram/matrix/"
outputLocG = "/home/egg/Projects/Papantonis_Stag/Results/differential_loop_o_gram/graphs/"

# Loop on hic files
for i in range(0,len(matrixFileNameList)):

  # Parameters
  matrixN = "3B9"
  if("127" in matrixFileNameList[i]): matrixN = "127"
  matrixSignal = "minus"
  if("plus" in matrixFileNameList[i]): matrixSignal = "plus"
  matrixName = ":".join(["matrix", matrixN, matrixSignal])
  regionN = regionFileNameList[i].split("/")[-1].split(".")[0].split("_")[1]
  regionType = regionFileNameList[i].split("/")[-1].split(".")[0].split("_")[2]
  regionName = ":".join(["region", regionN, regionType])
  finalName = matrixName+"_"+regionName

  # Input
  loopBins = "14"
  resolution = "10000"
  regionFileName = regionFileNameList[i]
  matrixFileName = matrixFileNameList[i]
  outputMatrixFileName = outputLocM+finalName+".txt"
  outputGraphFileName = outputLocG+finalName+".pdf"

  # Creating matrix
  command = "python loopOgram.py "+" ".join([loopBins, resolution, regionFileName, matrixFileName, outputMatrixFileName])
  os.system(command)

  # Creating heatmap
  command = "R CMD BATCH '--args '"+outputMatrixFileName+"' '"+outputGraphFileName+" loopOgramHeatmap.R loopOgramHeatmap.Rout"
  os.system(command)
    

