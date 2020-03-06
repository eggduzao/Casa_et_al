
# Import
import os
import sys

"""
# Stag List
counter = 1
fl = "/projects/ag-papan/eduardo/Wendt_Stag/Code/8_Meta_Plots/input/"
ml = "/projects/ag-papan/eduardo/Wendt_Stag/Data/stag_matrix_files/25K_norm/"
tl = "/projects/ag-papan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/"
ol = "/projects/ag-papan/eduardo/Wendt_Stag/Results/8_Meta_Plots/4_Meta_Tad_Tables/"
stagList = ["STAG1", "STAG2", "WT"]

# Stag Loop
for stag in stagList:

  # Diff Type List
  diffTypeList = ["0_equal", "1_split", "2_merge", "3_shift_upstream", "4_shift_downstream", "5_new"]

  # Diff Type Loop
  for diff in diffTypeList:

    # Signal List
    signalList = [["3B9_5-","."], ["3B9_5plus","."], ["69_127-","."], ["69_127plus","."],
                  ["69_127plus","69_127-"], ["3B9_5plus","3B9_5-"], ["69_127-","3B9_5-"], ["69_127plus","3B9_5plus"]]
    signalLabelList = [["STAG2-","."], ["STAG2+","."], ["STAG1-","."], ["STAG1+","."],
                  ["STAG1+","STAG1-"], ["STAG2+","STAG2-"], ["STAG1-","STAG2-"], ["STAG1+","STAG2+"]]

    # Signal Loop
    for i in range(0, len(signalList)):

      # Parameters
      signal = signalList[i]
      if("." in signalLabel): signalLabel = signalLabelList[i][0]
      else: signalLabel = "_vs_".join(signalLabelList[i])
      outLoc = ol + "/".join([stag, diff]) + "/"
      command = "mkdir -p "+outLoc
      os.system(command)

      # Input
      resolution = "25000"
      number_of_bins = "100"
      tad_1_file_name = tl + "_".join([diff, stag, "TAD1"]) + ".txt"
      tad_2_file_name = tl + "_".join([diff, stag, "TAD2"]) + ".txt"
      hic_1_file_name = ml + signal[0] + ".txt"
      hic_2_file_name = ml + signal[1] + ".txt"
      output_file_name = outLoc + signalLabel + ".txt"

      # Write
      inFile = open(fl + str(counter) + "_rmt.txt", "w")
      inFile.write("\n".join([resolution, number_of_bins, tad_1_file_name, tad_2_file_name, hic_1_file_name, hic_2_file_name, output_file_name]))
      inFile.close()
      counter += 1
"""

# /media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/5_new_STAG1_TAD1.txt
# /media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/5_new_STAG1_TAD2.txt

# Input
resolution = "25000"
number_of_bins = "10"
tad_1_file_name = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/5_new_STAG2_TAD1.txt"
tad_2_file_name = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Results/8_Meta_Plots/3_Tad_Files/5_new_STAG2_TAD1.txt"
hic_1_file_name = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/stag_matrix_files/25K_norm/3B9_5plus.txt"
hic_2_file_name = "/media/egg/cheops_agpapan/eduardo/Wendt_Stag/Data/stag_matrix_files/25K_norm/3B9_5-.txt"
output_file_name = "/home/egg/Desktop/metatad.txt"

command = "python 4_metaTAD_HicSignal.py "+" ".join([resolution, number_of_bins, tad_1_file_name, tad_2_file_name, hic_1_file_name, hic_2_file_name, output_file_name])
os.system(command)


