
# Import
import os
import sys

# Input List
counter = 1
ifn = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/29_loops/input/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/fithic_matrix_files/25K_norm/"
ol = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/29_loops/2_contact_files_raw/"
inList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]

# Input Loop
for inName in inList:

  # Chromosome List
  chromList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chromList:

    # lb List
    lbList = ["20000"]

    # lb Loop
    for lb in lbList:

      # ub List
      ubList = ["200000"]

      # ub Loop
      for ub in ubList:

        # oc List
        ocList = ["200"]

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
            outputLocation = ol + inName + "/"
            os.system("mkdir -p "+outputLocation)

            # Execution
            inFileName = ifn+str(counter)+".txt"
            inFile = open(inFileName, "w")
            inFile.write("\n".join([resolution, lowerBound, upperBound, equalOcBinNb, splineTimes, contactFileName, fragmentFileName, outputName, outputLocation]))
            inFile.close()
            counter += 1


