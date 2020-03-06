
# Import
import os
import sys

###################################################################################################
# Input
###################################################################################################

# Input
chromSizesFileName = sys.argv[1]
ensemblGeneFileName = sys.argv[2]
tempLocation = sys.argv[3]
outputFileName = sys.argv[4]

# Initialization
ext = 2000
command = "mkdir -p "+tempLocation
os.system(command)
  
###################################################################################################
# Creating Gene dictionaries
###################################################################################################

# Creating gene / exon / intron / tss / tes dictionaries
ensemblGeneFile = open(ensemblGeneFileName, "rU")
ensemblGeneFile.readline()
geneDict = dict()
tssDict = dict()
ttsDict = dict()
exonDict = dict()
intronDict = dict()
for line in ensemblGeneFile:

  ll = line.strip().split("\t")
  ensg = ll[12]
  score =ll[8]
  try:
    checkRepetition = geneDict[ensg]
    if(int(checkRepetition[4]) >= int(score)): continue
  except Exception: pass

  # Gene / tss / tes dictionary
  chrom = ll[2]; strand = ll[3]; tss = ll[4]; tes = ll[5]
  if(strand == "+"):
    start = ll[4]; end = ll[5]
    tss1 = str(max(int(ll[4]) - ext,0)); tss2 = str(int(ll[4]) + ext)
    tes1 = str(max(int(ll[5]) - ext,0)); tes2 = str(int(ll[5]) + ext)
  else:
    start = ll[4]; end = ll[5]
    tss1 = str(max(int(ll[5]) - ext,0)); tss2 = str(int(ll[5]) + ext)
    tes1 = str(max(int(ll[4]) - ext,0)); tes2 = str(int(ll[4]) + ext)
  try: geneSymbol = aliasDict[ensg]
  except Exception: geneSymbol = ensg
  geneDict[ensg] = [chrom, start, end, "GENE:"+geneSymbol, score, strand]
  tssDict[ensg] = [chrom, tss1, tss2, "TSS:"+geneSymbol, score, strand]
  ttsDict[ensg] = [chrom, tes1, tes2, "TTS:"+geneSymbol, score, strand]

  # Exon / Intron dictionaries
  eList1 = ll[9].split(",")[:-1]
  eList2 = ll[10].split(",")[:-1]
  if(int(eList1[0]) > int(tss)): intronDict[ensg] = [chrom, tss, eList1[0], "INTRON:"+geneSymbol, score, strand]
  prev = None
  for i in range(0,len(eList1)):
    if(prev): intronDict[ensg] = [chrom, prev, eList1[i], "INTRON:"+geneSymbol, score, strand]
    exonDict[ensg] = [chrom, eList1[i], eList2[i], "EXON:"+geneSymbol, score, strand]
    prev = eList2[i]
  if(int(eList2[-1]) < int(tes)): intronDict[ensg] = [chrom, tss, eList2[0], "INTRON:"+geneSymbol, score, strand]

ensemblGeneFile.close()

# Writing dictionaries
genomicTempFileName = tempLocation+"genomicTempFileName.bed"
genomicTempFile = open(genomicTempFileName, "w")
keys = geneDict.keys()
for k in keys:
  try: genomicTempFile.write("\t".join(geneDict[k])+"\n")
  except Exception: pass
  try: genomicTempFile.write("\t".join(tssDict[k])+"\n")
  except Exception: pass
  try: genomicTempFile.write("\t".join(ttsDict[k])+"\n")
  except Exception: pass
  try: genomicTempFile.write("\t".join(exonDict[k])+"\n")
  except Exception: pass
  try: genomicTempFile.write("\t".join(intronDict[k])+"\n")
  except Exception: pass
genomicTempFile.close()

###################################################################################################
# Intergenic Regions
###################################################################################################

# Reading chrom sizes
chromSizesDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1])
chromSizesFile.close()
chromList = chromSizesDict.keys()

# Iterating through each chromosome
for chrom in chromList:
  tempChromFileName = tempLocation+"tempChromFileName.bed"
  tempChromFile = open(tempChromFileName, "w")
  tempChromFile.write("\t".join([chrom, "0", str(chromSizesDict[chrom])])+"\n")
  tempChromFile.close()
  intersectTempFileName = tempLocation+"intersectTempFileName.bed"
  command = "subtractBed -a "+tempChromFileName+" -b "+genomicTempFileName+" > "+intersectTempFileName
  os.system(command)
  intersectTempFile = open(intersectTempFileName, "rU")
  genomicTempFile = open(genomicTempFileName, "a")
  for line in intersectTempFile:
    ll = line.strip().split("\t")
    genomicTempFile.write("\t".join(ll+["INTERGENIC:.", "0", "."])+"\n")
  intersectTempFile.close()
  genomicTempFile.close()

###################################################################################################
# Creating BAM file
###################################################################################################

# Grep chromosomes
grepFileName = tempLocation+"grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+genomicTempFileName+" > "+grepFileName
os.system(command)

# Sort bed
sortFileName = tempLocation+"sortFileName.bed"
command = "sort -k1,1 -k2,2n "+grepFileName+" > "+sortFileName
os.system(command)

# Bed To Bam
bamFileName = tempLocation+"bamFileName.bam"
command = "bedToBam -i "+sortFileName+" -g "+chromSizesFileName+" > "+bamFileName
os.system(command)

# Sort Bam
command = "samtools sort "+bamFileName+" -o "+outputFileName
os.system(command)

# Index Bam
command = "samtools index "+outputFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


