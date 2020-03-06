
# Import
import os
import sys

# Cond List
fl = "/usr/users/egadegu/Projects/Wendt_Stag/Code/36_Losada_TADs/input/"
ilf = "/usr/users/egadegu/Projects/Wendt_Stag/Results/26_Losada_Interaction_Decay/1_Sparse_Matrices/25K_norm/"
olf = "/usr/users/egadegu/Projects/Wendt_Stag/Results/36_Losada_TADs/1_TADs_GMAP/"
condList = ["LCONT", "SISA1", "SISA2"]
res = "25000"

# Opening file
inFileName = fl + "1_gmap.txt"
inFile = open(inFileName,"w")

# Input Loop
for cond in condList:

  # Chromosome List
  chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

  # Chromosome Loop
  for chrom in chrList:

    # Logt List
    logtList = ["T"]

    # Logt Loop
    for logt in logtList:

      # DomOrder List
      domorderList = ["2"]

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
                  hthrList = ["0.99"]

                  # Hthr Loop
                  for hthr in hthrList:

                    # T1thr List
                    t1thrList = ["0.75"]

                    # T1thr Loop
                    for t1thr in t1thrList:
              
                      # Name
                      namy = cond + "/" + chrom
                      name = cond + "_" + chrom

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
                      inputFileName = ilf + namy + "_bin.txt"
                      outputFilePrefix = olf + name

                      # Execution
                      #inFile.write("'--args '" + "' '".join([resolution, logt, domOrder, maxDistInBin, minD, maxD, minDP, maxDP, hthr, t1thr, inputFileName, outputFilePrefix])+"\n")
                      inFile.write(" ".join([resolution, logt, domOrder, maxDistInBin, minD, maxD, minDP, maxDP, hthr, t1thr, inputFileName, outputFilePrefix])+"\n")

inFile.close()


