
# Import
import os
import sys

# Input List
counter = 1
ifn = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/5_Loops/input/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/fithic_matrix_files/25K_norm/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/5_Loops/2_contact_files_raw/"
inList = ["69_127-", "69_127plus", "3B9_5-", "3B9_5plus"]
labelList = ["STAG1_minus_auxin", "STAG1_plus_auxin", "STAG2_minus_auxin", "STAG2_plus_auxin"]

# Input Loop
for i in range(0,len(inList)):

  # Names
  inName = inList[i]
  inLabel = labelList[i]

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # lb List
    lbList = ["50000", "150000"]

    # lb Loop
    for lb in lbList:

      # ub List
      ubList = ["5000000", "10000000"]

      # ub Loop
      for ub in ubList:

        # oc List
        ocList = ["100"]

        # oc Loop
        for oc in ocList:

          # sp List
          spList = ["1"]

          # sp Loop
          for sp in spList:

            # Input
            resolution = "25000"
            lowerBound = lb
            upperBound = ub
            equalOcBinNb = oc
            splineTimes = sp
            contactFileName = il + "_".join([inName, chrom, "contact.gz"])
            fragmentFileName = il + "_".join([inName, chrom, "fragment.gz"])
            outputName = chrom
            outputLocation = ol + "_".join([inLabel, lb, ub]) + "/"
            os.system("mkdir -p "+outputLocation)

            # Execution
            inFileName = ifn + str(counter) + "_rfh.txt"
            inFile = open(inFileName, "w")
            inFile.write("\n".join([resolution, lowerBound, upperBound, equalOcBinNb, splineTimes, contactFileName, fragmentFileName, outputName, outputLocation]))
            inFile.close()
            counter += 1


