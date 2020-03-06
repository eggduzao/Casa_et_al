
# Import
import os
import sys

# Cond List
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Code/24_TADs_multiple_parameter_test/input/"
ilf = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/stag_matrix_files/25K_norm/"
olf = "/projects/ag-papan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/"
condList = ["3B9_5-", "3B9_5plus", "69_127-", "69_127plus"]
res = "25000"
counter = 1

# Input Loop
for cond in condList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Logt List
    logtList = ["T", "F"]

    # Logt Loop
    for logt in logtList:

      # DomOrder List
      domorderList = ["2", "3"]

      # DomOrder Loop
      for domorder in domorderList:

        # MaxDist List
        maxDistList = ["80"]

        # MaxDist Loop
        for maxDist in maxDistList:

          # Mind List
          mindList = ["25"]

          # Mind Loop
          for mind in mindList:

            # Maxd List
            maxdList = ["100"]

            # Maxd Loop
            for maxd in maxdList:

              # Mindp List
              mindpList = ["5"]

              # Mindp Loop
              for mindp in mindpList:

                # Maxdp List
                maxdpList = ["10"]

                # Maxdp Loop
                for maxdp in maxdpList:

                  # Hthr List
                  hthrList = ["0.95", "0.99"]

                  # Hthr Loop
                  for hthr in hthrList:

                    # T1thr List
                    t1thrList = ["0.5", "0.75"]

                    # T1thr Loop
                    for t1thr in t1thrList:
              
                      # Name
                      namy = cond+"/"+chrom
                      outLoc = olf+namy+"/"
                      os.system("mkdir -p "+outLoc)
                      name = "_".join([logt, domorder, maxDist, mind, maxd, mindp, maxdp, hthr, t1thr])

                      # Input
                      resolution = res
                      logt = logt
                      domOrder = domorder
                      maxDistInBin = maxDist
                      minD = mind
                      maxD = maxd
                      minDP = mindp
                      maxDP = maxdp
                      hthr = hthr
                      t1thr = t1thr
                      inputFileName = ilf+namy+"_bin.txt"
                      outputFilePrefix = outLoc+name

                      # Execution
                      inFileName = il+str(counter)+".txt"
                      inFile = open(inFileName, "w")
                      inFile.write("\n".join([resolution, logt, domOrder, maxDistInBin, minD, maxD, minDP, maxDP, hthr, t1thr, inputFileName, outputFilePrefix]))
                      inFile.close()
                      counter += 1


