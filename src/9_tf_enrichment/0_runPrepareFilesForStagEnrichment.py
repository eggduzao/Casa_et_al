
# Import
import os
import sys

# Peak List
il1 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/macs/"
il2 = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/footprints/"
ol = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/tf_enrichment_new/input/"
peakList = [il1+"DNase_peaks_filter.narrowPeak", il1+"DNase_summits_filter.bed", il1+"H3K27ac_peaks_filter.narrowPeak", il1+"H3K27ac_summits_filter.bed", il2+"footprints.bed"]

# Peak Loop
for peakFile in peakList:

  # Stag List
  sl = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Data/stag_bed_files/"
  stagList = [sl+"STAG1only_mrg_replicates.bed", sl+"STAG2only_mtg_replicates.bed"]

  # Stag Loop
  for stagFile in stagList:

    # Parameters
    pname = peakFile.split("/")[-1].split(".")[0]
    if("STAG1" in stagFile.split("/")[-1].split(".")[0]): sname = "stag1"
    else: sname = "stag2"
    name = "_".join([pname, sname])
    fext = peakFile.split("/")[-1].split(".")[-1]

    # Input
    if(fext == "narrowPeak"): ext = "500"
    elif("summits" in pname): ext = "500"
    else: ext = "25"
    peakFileName = peakFile
    stagFileName = stagFile
    tempLocation = ol+"TEMP/"
    outputFileName = ol+name+".bed"

    # Execution
    command = "python prepareFilesForStagEnrichment.py "+" ".join([ext, peakFileName, stagFileName, tempLocation, outputFileName])
    os.system(command)


