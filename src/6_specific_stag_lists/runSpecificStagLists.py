
# Import
import os
import sys

# Stag List
il = "/home/egg/Projects/Papantonis_Stag/Results/ranked_stag_lists/"
ol = "/home/egg/Projects/Papantonis_Stag/Results/specific_stag_lists/"
stagList = [["69_127-_STAG1only/69_127-_STAG1only_full.txt", "69_127plus_STAG1only/69_127plus_STAG1only_full.txt"],
            ["3B9_5-_STAG2only/3B9_5-_STAG2only_full.txt", "3B9_5plus_STAG2only/3B9_5plus_STAG2only_full.txt"]]
stagLabelList = ["STAG1", "STAG2"]

# Stag Loop
for i in range(0,len(stagList)):

  # Input
  plusAuxinFileName = il+stagList[i][1]
  minusAuxinFileName = il+stagList[i][0]
  specificFileName = ol+"auxin_specific/"+stagLabelList[i]+".txt"
  sensitiveFileName = ol+"auxin_sensitive/"+stagLabelList[i]+".txt"
  bothFileName = ol+"both/"+stagLabelList[i]+".txt"

  # Execution
  command = "python specificStagLists.py "+" ".join([plusAuxinFileName, minusAuxinFileName, specificFileName, sensitiveFileName, bothFileName])
  os.system(command)
    
