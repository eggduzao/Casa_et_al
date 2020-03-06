
# Import
import os
import sys

# Input location
il = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/1_TADs_GMAP/"
outTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/0_total_TAD_counts/number_tads.txt"
outhTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/0_total_TAD_counts/number_htads.txt"
outthTadFileName = "/media/egg/cheops_agpapan/eduardo/Papantonis_Stag/Results/24_TADs_multiple_parameter_test/0_total_TAD_counts/number_all.txt"

# Parameters
header = ["CASE", "STAG1+", "STAG1-", "STAG2+", "STAG2-"]
outTadFile = open(outTadFileName,"w")
outhTadFile = open(outhTadFileName,"w")
outthTadFile = open(outthTadFileName,"w")
outTadFile.write("\t".join(header)+"\n")
outhTadFile.write("\t".join(header)+"\n")
outthTadFile.write("\t".join(header)+"\n")

# Case List
caseList = ["F_2_80_25_100_5_10_0.95_0.5", "F_2_80_25_100_5_10_0.95_0.75", "F_2_80_25_100_5_10_0.99_0.5", "F_2_80_25_100_5_10_0.99_0.75",        "F_3_80_25_100_5_10_0.95_0.5", "F_3_80_25_100_5_10_0.95_0.75", "F_3_80_25_100_5_10_0.99_0.5", "F_3_80_25_100_5_10_0.99_0.75", "T_2_80_25_100_5_10_0.95_0.5", "T_2_80_25_100_5_10_0.95_0.75", "T_2_80_25_100_5_10_0.99_0.5", "T_2_80_25_100_5_10_0.99_0.75", "T_3_80_25_100_5_10_0.95_0.5", "T_3_80_25_100_5_10_0.95_0.75", "T_3_80_25_100_5_10_0.99_0.5", "T_3_80_25_100_5_10_0.99_0.75"]
caseNumberList = [str(e) for e in range(1,len(caseList)+1)]

# Case Loop
for k in range(0,len(caseList)):

  case = caseList[k]
  caseN = caseNumberList[k]

  # Name List
  nameList = ["69_127plus", "69_127-", "3B9_5plus", "3B9_5-"]
  ttotalCount = [0, 0, 0, 0]
  htotalCount = [0, 0, 0, 0]
  thtotalCount = [0, 0, 0, 0]

  # Chromosome Loop
  for chrom in ["chr"+str(e) for e in range(1,23)+["X"]]:

    for i in range(0,len(nameList)):
      inFileName = il+nameList[i]+"/"+chrom+"/"+case+"_htad.txt"
      inFile = open(inFileName,"rU")
      inFile.readline()
      tadCounter = 0
      htadCounter = 0
      for line in inFile:
        ll = line.strip().split("\t")
        if(ll[2] == "1"):
          ttotalCount[i] += 1
          thtotalCount[i] += 1
        else:
          htotalCount[i] += 1
          thtotalCount[i] += 1
      inFile.close()

  outTadFile.write("\t".join([caseN]+[str(e) for e in ttotalCount])+"\n")
  outhTadFile.write("\t".join([caseN]+[str(e) for e in htotalCount])+"\n")
  outthTadFile.write("\t".join([caseN]+[str(e) for e in thtotalCount])+"\n")

# Closing files
outTadFile.close()
outhTadFile.close()
outthTadFile.close()


