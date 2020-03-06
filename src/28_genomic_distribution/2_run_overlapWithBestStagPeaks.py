
# Import
import os
import sys

# Stag List
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
sl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/hct116_stag_signal/"
stagList = ["STAG1", "STAG2"]
stagCount = ["370424512", "355785002"]

# Stag Loop
for i in range(0,len(stagList)):

  st = stagList[i]
  sc = stagCount[i]

  # Input
  stagTotalCount = sc
  ctcfTotalCount = "25415806"
  stagFileName = il+st+"bestpeaks_filter.bed"
  if(st == "STAG1"): stagBwSignalFileName = sl+st+"_127_EGFP.bw"
  else: stagBwSignalFileName = sl+st+"_3B9_EGFP.bw"
  ctcfSignalFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/tf/bam/HCT116_ChIP-seq_CTCF_BROAD.bam"
  genomicRegionsFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/1_genomic_regions/regions.bam"
  chromHmmRegionsFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/16_chromhmm/HCT116_20_dense.bam"
  outputFileNamePrefix = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/28_genomic_distribution/2_stag_regions_overlap/"+st

  # Execution
  command = "python 2_overlapWithBestStagPeaks.py "+" ".join([stagTotalCount, ctcfTotalCount, stagFileName, stagBwSignalFileName, ctcfSignalFileName, genomicRegionsFileName, chromHmmRegionsFileName, outputFileNamePrefix])
  os.system(command)


