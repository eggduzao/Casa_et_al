
# Import
import os
import sys
from pysam import Samfile

###################################################################################################
# Input
###################################################################################################

# Input
aliasFileName = sys.argv[1]
chromSizesFileName = sys.argv[2]
ensemblDictFileName = sys.argv[3]
ensemblGeneFileName = sys.argv[4]
ensemblBamFileName = sys.argv[5]
ensemblNoIEBedFileName = sys.argv[6]
ensemblNoIEBamFileName = sys.argv[7]
enhancerFileName = sys.argv[8]
h3k4me1FileName = sys.argv[9]
h3k4me3FileName = sys.argv[10]
h3k27acFileName = sys.argv[11]
ctcfFileName = sys.argv[12]
tempLocation = sys.argv[13]
outputBedFileName = sys.argv[14]
outputBamFileName = sys.argv[15]

# Initialization
ext = 2000
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Creating chromsizes and alias dictionaries
###################################################################################################

# Alias dictionary
aliasDict = dict()
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g.upper()] = value.upper()
aliasFile.close()

# Reading chrom sizes
chromSizesDict = dict()
chromSizesFile = open(chromSizesFileName,"rU")
for line in chromSizesFile:
  ll = line.strip().split("\t")
  chromSizesDict[ll[0]] = int(ll[1])
chromSizesFile.close()
chromList = chromSizesDict.keys()

# Creating ENSEMBL dictionary file
ensemblDict = dict() # ENST -> ENSG
ensemblDictFile = open(ensemblDictFileName, "rU")
for line in ensemblDictFile:
  ll = line.strip().split("\t")
  ensemblDict[ll[7]] = ll[0]
ensemblDictFile.close()

###################################################################################################
# Creating Gene dictionaries
###################################################################################################

# Initialization
h3k4me3File = Samfile(h3k4me3FileName, "rb")
h3k27acFile = Samfile(h3k27acFileName, "rb")

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
  if(chrom not in chromList): continue
  if(strand == "+"):
    start = ll[4]; end = ll[5]
    tss1 = str(max(int(ll[4]) - ext,0)); tss2 = str(int(ll[4]) + 0)
    tes1 = str(max(int(ll[5]) - 0,0)); tes2 = str(int(ll[5]) + ext)
    p1 = tss1; p2 = end
  else:
    start = ll[4]; end = ll[5]
    tss1 = str(max(int(ll[5]) - 0,0)); tss2 = str(int(ll[5]) + ext)
    tes1 = str(max(int(ll[4]) - ext,0)); tes2 = str(int(ll[4]) + 0)
    p1 = start; p2 = tss2
  try: geneSymbol = aliasDict[ensg]
  except Exception: geneSymbol = ensg

  # Active status
  activeStatus = "INACTIVE"
  h3k4me3Fetch = h3k4me3File.fetch(chrom, int(p1), int(p2))
  mSum = sum(1 for _ in h3k4me3Fetch)
  h3k27acFetch = h3k27acFile.fetch(chrom, int(p1), int(p2))
  aSum = sum(1 for _ in h3k27acFetch)
  if(mSum > 0 or aSum > 0): activeStatus = "ACTIVE"

  geneDict[ensg] = [chrom, start, end, "GENE:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
  tssDict[ensg] = [chrom, tss1, tss2, "PROMOTER:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
  ttsDict[ensg] = [chrom, tes1, tes2, "TTS:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]

  # Exon / Intron dictionaries
  eList1 = ll[9].split(",")[:-1]
  eList2 = ll[10].split(",")[:-1]
  if(int(eList1[0]) > int(tss)): intronDict[ensg] = [chrom, tss, eList1[0], "INTRON:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
  prev = None
  for i in range(0,len(eList1)):
    if(prev): intronDict[ensg] = [chrom, prev, eList1[i], "INTRON:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
    exonDict[ensg] = [chrom, eList1[i], eList2[i], "EXON:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
    prev = eList2[i]
  if(int(eList2[-1]) < int(tes)): intronDict[ensg] = [chrom, tss, eList2[0], "INTRON:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]

ensemblGeneFile.close()

# Genes without TSS/EXON
ensemblNoIEBedFile = open(ensemblNoIEBedFileName, "rU")
for line in ensemblNoIEBedFile:

  ll = line.strip().split("\t")
  chrom = ll[0]; start = ll[1]; end = ll[2]; nn = ll[3].split(":"); score = ll[4]; strand = ll[5]
  ensg = nn[0]; name = nn[1]; gene = nn[2]
  if(gene == "."): gene = name

  try:
    checkRepetition = geneDict[ensg]
    continue
  except Exception: pass

  # Gene / tss / tes dictionary
  if(chrom not in chromList): continue
  if(strand == "+"):
    gene_start = start; gene_end = end
    tss1 = str(max(int(gene_start) - ext, 0)); tss2 = str(int(gene_start) + 0)
    tes1 = str(max(int(gene_end) - 0, 0)); tes2 = str(int(gene_end) + ext)
    p1 = tss1; p2 = gene_end
  else:
    gene_start = start; gene_end = end
    tss1 = str(max(int(gene_end) - 0,0)); tss2 = str(int(gene_end) + ext)
    tes1 = str(max(int(gene_start) - ext,0)); tes2 = str(int(gene_start) + 0)
    p1 = gene_start; p2 = tss2
  try: geneSymbol = aliasDict[ensg]
  except Exception:
    if(gene == "."): geneSymbol = gene
    geneSymbol = ensg

  # Active status
  activeStatus = "INACTIVE"
  h3k4me3Fetch = h3k4me3File.fetch(chrom, int(p1), int(p2))
  mSum = sum(1 for _ in h3k4me3Fetch)
  h3k27acFetch = h3k27acFile.fetch(chrom, int(p1), int(p2))
  aSum = sum(1 for _ in h3k27acFetch)
  if(mSum > 0 or aSum > 0): activeStatus = "ACTIVE"

  geneDict[ensg] = [chrom, start, end, "GENE:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
  tssDict[ensg] = [chrom, tss1, tss2, "PROMOTER:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]
  ttsDict[ensg] = [chrom, tes1, tes2, "TTS:"+geneSymbol+":"+activeStatus, str(int(float(score))), strand]

ensemblNoIEBedFile.close()

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

###################################################################################################
# Enhancer Regions with EnhancerAtlas
###################################################################################################

# Initialization
step = 1000000
enhancersFile = Samfile(enhancerFileName, "rb")

# Iterating through each chromosome
for chrom in chromList:

  # Iterating through genome
  for i in range(0, chromSizesDict[chrom], step):

    # Region
    p1 = i; p2 = min(chromSizesDict[chrom], i + step)
    if(p1 >= chromSizesDict[chrom]-1): continue

    # Iterating through enhancers
    enhancerFetch = enhancersFile.fetch(chrom, p1, p2)
    for read in enhancerFetch:

      # Read coordintes
      startx = str(read.pos)
      endx = str(read.aend)

      # Active status
      activeStatus = "POISED"
      h3k27acFetch = h3k27acFile.fetch(chrom, int(startx), int(endx))
      aSum = sum(1 for _ in h3k27acFetch)
      if(aSum > 0): activeStatus = "ACTIVE"

      # Name and score
      nameVec = read.qname.split(":")
      score = nameVec[2]
      try: gene = ensemblDict[nameVec[1]]
      except Exception: gene = nameVec[1]
      try: geneSymbol = aliasDict[gene]
      except Exception: geneSymbol = gene

      # Writing enhancer
      genomicTempFile.write("\t".join([chrom, startx, endx, "ENHANCER:"+geneSymbol+":"+activeStatus, str(int(float(score))), "+"])+"\n")

###################################################################################################
# Enhancer Regions with H3K4me1+H3K27ac
###################################################################################################

# Initialization
ensemblBamFile = Samfile(ensemblBamFileName, "rb")
ensemblBamFile2 = Samfile(ensemblNoIEBamFileName, "rb")
h3k4me1File = Samfile(h3k4me1FileName, "rb")

# Iterating through each chromosome
for chrom in chromList:

  # Iterating through genome
  for i in range(0, chromSizesDict[chrom], step):

    # Region
    p1 = i; p2 = min(chromSizesDict[chrom], i + step)
    if(p1 >= chromSizesDict[chrom]-1): continue

    # Iterating on H3K4me1
    h3k4me1Fetch = h3k4me1File.fetch(chrom, p1, p2)
    for read in h3k4me1Fetch:

      # Read coordintes
      startx = str(read.pos)
      endx = str(read.aend)

      # Check if enhancer or if it already exists
      enhancerFetch = enhancersFile.fetch(chrom, int(startx), int(endx))
      eSum = sum(1 for _ in enhancerFetch)
      ensemblFetch = ensemblBamFile.fetch(chrom, int(startx), int(endx))
      gSum1 = sum(1 for _ in ensemblFetch)
      ensemblFetch2 = ensemblBamFile2.fetch(chrom, int(startx), int(endx))
      gSum2 = sum(1 for _ in ensemblFetch2)
      gSum = gSum1 + gSum2
      h3k4me3Fetch = h3k4me3File.fetch(chrom, int(startx), int(endx))
      m3Sum = sum(1 for _ in h3k4me3Fetch)
      if(eSum > 0 or gSum > 0 or m3Sum > 0): continue

      # Active status
      activeStatus = "POISED"
      h3k27acFetch = h3k27acFile.fetch(chrom, int(startx), int(endx))
      aSum = sum(1 for _ in h3k27acFetch)
      if(aSum > 0): activeStatus = "ACTIVE"   

      # Writing enhancer
      genomicTempFile.write("\t".join([chrom, startx, endx, "ENHANCER:.:"+activeStatus, "1000", "+"])+"\n")
 
    # Iterating on H3K27ac
    h3k27acFetch = h3k27acFile.fetch(chrom, p1, p2)
    for read in h3k4me1Fetch:

      # Read coordintes
      startx = str(read.pos)
      endx = str(read.aend)

      # Check if enhancer or if it already exists
      enhancerFetch = enhancersFile.fetch(chrom, int(startx), int(endx))
      eSum = sum(1 for _ in enhancerFetch)
      ensemblFetch = ensemblBamFile.fetch(chrom, int(startx), int(endx))
      gSum1 = sum(1 for _ in ensemblFetch)
      ensemblFetch2 = ensemblBamFile2.fetch(chrom, int(startx), int(endx))
      gSum2 = sum(1 for _ in ensemblFetch2)
      gSum = gSum1 + gSum2
      h3k4me3Fetch = h3k4me3File.fetch(chrom, int(startx), int(endx))
      m3Sum = sum(1 for _ in h3k4me3Fetch)
      if(eSum > 0 or gSum > 0 or m3Sum > 0): continue  

      # Writing enhancer
      genomicTempFile.write("\t".join([chrom, startx, endx, "ENHANCER:.:ACTIVE", "1000", "+"])+"\n")
      
# Closing files
genomicTempFile.close()
ensemblBamFile.close()
enhancersFile.close()
h3k4me1File.close()
h3k4me3File.close()
h3k27acFile.close()

###################################################################################################
# Superenhancer Regions
###################################################################################################

# Initialization
maxDist = "12500"

# Fetching enhancers
enhancerTempFileName = tempLocation+"enhancerTempFileName.bed"
command = "grep ENHANCER "+genomicTempFileName+" > "+enhancerTempFileName
os.system(command)

# Uniq enhancers
uniqEnhTempFileName = tempLocation+"uniqEnhTempFileName.bed"
command = "sort "+enhancerTempFileName+" | uniq > "+uniqEnhTempFileName
os.system(command)

# Sorting enhancers
enhancerTempSortFileName = tempLocation+"enhancerTempSortFileName.bed"
command = "sort -k1,1 -k2,2n "+uniqEnhTempFileName+" > "+enhancerTempSortFileName
os.system(command)

# Fetching superenhancers
superEnhFileName = tempLocation+"superEnhFileName.bed"
command = "mergeBed -d "+maxDist+" -c 4,5,6 -o collapse,mean,distinct -i "+enhancerTempSortFileName+" > "+superEnhFileName
os.system(command)

# Writing superenhancers
superEnhFile = open(superEnhFileName, "rU")
genomicTempFile = open(genomicTempFileName, "a")
for line in superEnhFile:
  ll = line.strip().split("\t")
  name = ll[3].split(",")
  if(len(name) <= 3): continue
  actStatus = name[0].split(":")[2]
  genomicTempFile.write("\t".join([ll[0], ll[1], ll[2], "SUPERENHANCER:.:"+actStatus, "1000", "+"])+"\n")

# Closing files
superEnhFile.close()

###################################################################################################
# CTCF regions
###################################################################################################

# Initialization
step = 1000000
ctcfFile = Samfile(ctcfFileName, "rb")

# Iterating through each chromosome
for chrom in chromList:

  # Iterating through genome
  for i in range(0, chromSizesDict[chrom], step):

    # Region
    p1 = i; p2 = min(chromSizesDict[chrom], i + step)
    if(p1 >= chromSizesDict[chrom]-1): continue

    # Iterating through enhancers
    ctcfFetch = ctcfFile.fetch(chrom, p1, p2)
    for read in ctcfFetch:

      # Read coordintes
      startx = str(read.pos)
      endx = str(read.aend)

      # Writing enhancer
      genomicTempFile.write("\t".join([chrom, startx, endx, "CTCF:.:.", "1000", "+"])+"\n")

# Closing files
ctcfFile.close()
genomicTempFile.close()

###################################################################################################
# Intergenic Regions
###################################################################################################

# Sorting genomicTempFile
genomicSortFileName = tempLocation+"genomicSortFileName.bed"
command = "sort -k1,1 -k2,2n "+genomicTempFileName+" > "+genomicSortFileName
os.system(command)

# Iterating through each chromosome
for chrom in chromList:
  tempChromFileName = tempLocation+"tempChromFileName.bed"
  tempChromFile = open(tempChromFileName, "w")
  tempChromFile.write("\t".join([chrom, "0", str(chromSizesDict[chrom])])+"\n")
  tempChromFile.close()
  intersectTempFileName = tempLocation+"intersectTempFileName.bed"
  command = "subtractBed -a "+tempChromFileName+" -b "+genomicSortFileName+" > "+intersectTempFileName
  os.system(command)
  intersectTempFile = open(intersectTempFileName, "rU")
  genomicTempFile = open(genomicSortFileName, "a")
  for line in intersectTempFile:
    ll = line.strip().split("\t")
    genomicTempFile.write("\t".join(ll+["INTERGENIC:.:.", "0", "+"])+"\n")
  intersectTempFile.close()
  genomicTempFile.close()

###################################################################################################
# Creating BAM file
###################################################################################################

# Grep chromosomes
grepFileName = tempLocation+"grepFileName.bed"
command = "grep -E 'chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX' "+genomicSortFileName+" > "+grepFileName
os.system(command)

# Uniq bed
uniqFileName = tempLocation+"uniqFileName.bed"
command = "sort "+grepFileName+" | uniq > "+uniqFileName
os.system(command)

# Sort bed
command = "sort -k1,1 -k2,2n "+uniqFileName+" > "+outputBedFileName
os.system(command)

# Bed To Bam
bamFileName = tempLocation+"bamFileName.bam"
command = "bedToBam -i "+outputBedFileName+" -g "+chromSizesFileName+" > "+bamFileName
os.system(command)

# Sort Bam
command = "samtools sort "+bamFileName+" -o "+outputBamFileName
os.system(command)

# Index Bam
command = "samtools index "+outputBamFileName
os.system(command)

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


