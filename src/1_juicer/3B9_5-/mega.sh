#!/bin/bash
# Makes a mega map from the files given on the command line, presumes hg19 
# MboI
# Note: this is an older version and should be updated

tmpDir="/projects/ag-papan/eduardo/Papantonis_Stag/Code/1_juicer/3B9_5-/TEMP"
juiceDir="/projects/ag-papan/eduardo/Papantonis_Stag/Code/1_juicer/juicerFolder"
genomeID="hg19"
site="NcoI"

# check usage
if [ $# -lt 2 ]
    then
    echo "Usage: `basename $0` <file1> <file2> ... <fileN>"
    echo " Combines the files on command line into mega map"
    exit 1
fi

# check that files exist, add to string
str=""
for var in "$@"
do
  if [ ! -f $var ]
      then
      echo "File $var does not exist"
      exit 1
  else
      str+="$var "
  fi
done

sort -T ${tmpDir} -k2,2d -k6,6d -m ${str} > merged_nodups.txt
if [ $? -eq 0 ]; then
    echo "merged_nodups successfully created"
else
    echo "merged_nodups failed, exiting"
    exit 1
fi
${juiceDir}/scripts/common/statistics.pl -q 1 -ointer.txt -s ${juiceDir}/restriction_sites/${genomeID}_${site}.txt merged_nodups.txt
if [ $? -eq 0 ]; then
    echo "inter.txt successfully created"
else
    echo "inter.txt failed, exiting"
    exit 1
fi
${juiceDir}/scripts/common/statistics.pl -q 30 -ointer_30.txt -s ${juiceDir}/restriction_sites/${genomeID}_${site}.txt merged_nodups.txt
if [ $? -eq 0 ]; then
    echo "inter_30.txt successfully created"
else
    echo "inter_30.txt failed, exiting"
    exit 1
fi
${juiceDir}/scripts/common/juicer_tools pre -f ${juiceDir}/restriction_sites/${genomeID}_${site}.txt -s inter.txt -g inter_hists.m -q 1 merged_nodups.txt inter.hic ${genomeID}
if [ $? -eq 0 ]; then
    echo "inter.hic successfully created"
else
    echo "inter.hic failed, exiting"
    exit 1
fi
${juiceDir}/scripts/common/juicer_tools pre -f ${juiceDir}/restriction_sites/${genomeID}_${site}.txt -s inter_30.txt -q 30 -g inter_30_hists.m merged_nodups.txt inter_30.hic ${genomeID}
if [ $? -eq 0 ]; then
    echo "inter_30.txt successfully created"
else
    echo "inter_30.txt failed, exiting"
    exit 1
fi


