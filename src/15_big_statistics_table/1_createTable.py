
# Import
import os
import sys
from pysam import Samfile, Fastafile

###################################################################################################
# Input
###################################################################################################

# Input
fileIsLoop = sys.argv[1]
stagFileName = sys.argv[2]
aliasFileName = sys.argv[3]
chromSizesFileName = sys.argv[4]
genomeFileName = sys.argv[5]
regionsFileName = sys.argv[6]
ensemblDictFileName = sys.argv[7]
enhancersFileName = sys.argv[8]
chrommHmmFileName = sys.argv[9]
expressionLabelList = sys.argv[10].split(",")
expressionFileNameList = sys.argv[11].split(",")
signalLabelList = sys.argv[12].split(",")
signalCountList = [int(e) for e in sys.argv[13].split(",")]
signalFileNameList = sys.argv[14].split(",")
controlCountList = [int(e) for e in sys.argv[15].split(",")]
controlFileNameList = sys.argv[16].split(",")
peakFileNameList = sys.argv[17].split(",")
motifLabelList = sys.argv[18].split(",")
motifFileNameList = sys.argv[19].split(",")
tempLocation = sys.argv[20]
outputFileName1 = sys.argv[21]
outputFileName2 = sys.argv[22]
outputFileName3 = sys.argv[23]
outputFileName4 = sys.argv[24]

# Initialization
geneExt = 5000
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def fileLen(fileName):
  if(os.stat(fileName).st_size <= 0): return 0
  i = 0
  with open(fileName) as f:
    for i, l in enumerate(f): pass
  return i + 1

def fetchTotalReadsBam(bamFile, region):
  returnN = 0
  for read in bamFile.fetch(region[0], region[1], region[2]): returnN += 1
  return returnN

def fetchSignalBam(bamFile, region, ext):
  regionLen = (region[2] - region[1])
  returnVec = [0.0] * regionLen
  for read in bamFile.fetch(region[0], region[1], region[2]):
    vecReadLoc = read.pos - region[1]
    for i in range(max(vecReadLoc - ext, 0), min(vecReadLoc + ext, len(returnVec))): returnVec[i] += 1.0
  return returnVec
  
###################################################################################################
# Creating Gene dictionaries
###################################################################################################

# Valid chromosomes
chrList = ["chr"+str(e) for e in range(1,23)+["X"]]

# Reading alias
aliasDict = dict() # alias -> gene_symbol
aliasFile = open(aliasFileName,"rU")
for line in aliasFile:
  ll = line.strip().split("\t")
  value = ll[1]
  geneList = [ll[0],ll[1]]+ll[2].split("&")
  for g in geneList: aliasDict[g] = value
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

# Creating GO dictionary file
goDict = dict() # ENSG -> GO
ensemblDictFile = open(ensemblDictFileName, "rU")
for line in ensemblDictFile:
  ll = line.strip().split("\t")
  goDict[ll[0]] = ll[4]
ensemblDictFile.close()

# Creating expression dictionaries
exprDictList = [dict() for e in range(0,len(expressionLabelList))] # ENSG -> [log2FC, minusAux, plusAux]
for i in range(0,len(expressionLabelList)):
  exprLabel = expressionLabelList[i]
  exprFileName = expressionFileNameList[i]
  exprFile = open(exprFileName, "rU")
  exprFile.readline()
  for line in exprFile:
    ll = line.strip().split(",")
    ensg = ll[0].replace("\"","")
    log2FoldChange = ll[3]
    minusAux = str(round((float(ll[8]) + float(ll[9]))/2,2))
    plusAux = str(round((float(ll[10]) + float(ll[11]))/2,2))
    exprDictList[i][ensg] = [log2FoldChange, minusAux, plusAux]
  exprFile.close()

###################################################################################################
# Execution
###################################################################################################

# Open bam files
genomeFile = Fastafile(genomeFileName)
regionsFile = Samfile(regionsFileName, "rb")
#chrommHmmFile = Samfile(chrommHmmFileName, "rb")
enhancersFile = Samfile(enhancersFileName, "rb")
signalFileList = [Samfile(e, "rb") for e in signalFileNameList]
controlFileList = [Samfile(e, "rb") for e in controlFileNameList]
motifFileList = [Samfile(e, "rb") for e in motifFileNameList]

# Creating RPM list
rpmList = [1000000./e for e in signalCountList]
rpmControlList = [1000000./e for e in controlCountList]

# Fetching index of CTCF motifs
ctcfIndexList = []
for i in range(0,len(motifLabelList)):
  if("CTCF_" in motifLabelList[i]): ctcfIndexList.append(i)

# Creating file header
stagFile = open(stagFileName, "rU")
headerRep = 1
if(fileIsLoop == "Y"): headerRep = 2
header1 = ["STAG_PEAK_CHROMOSOME", "STAG_PEAK_START", "STAG_PEAK_END", "STAG_PEAK_LENGTH", "STAG_PEAK_CG_CONTENT", "CTCF_MOTIF", "GENOMIC_REGION", "CHROMHMM_STATE", "CLOSEST_GENES_10K", "DISTANCE_TO_CLOSEST_GENES_TSS"] + [e+"_CLOSEST_GENES_EXPRESSION_MINUSAUX_PLUSAUX" for e in expressionLabelList] + ["CLOSEST_ENHANCERS_10K", "DISTANCE_TO_CLOSEST_ENHANCERS_CENTER"] + [e+"_CLOSEST_ENHANCERS_EXPRESSION_MINUSAUX_PLUSAUX" for e in expressionLabelList] + ["CLOSEST_GENES_GO_TERMS"]
header2 = ["STAG_PEAK_CHROMOSOME", "STAG_PEAK_START", "STAG_PEAK_END", "STAG_PEAK_LENGTH"] + ["TOTAL_"+e+"_NORMALIZED_COUNT(RPM)\t"+"AVERAGE_"+e+"_NORMALIZED_COUNT(RPKM)\tTOTAL_CONTROL_"+e+"_NORMALIZED_COUNT(RPM)\t"+"AVERAGE_CONTROL_"+e+"_NORMALIZED_COUNT(RPKM)" for e in signalLabelList]
header3 = ["STAG_PEAK_CHROMOSOME", "STAG_PEAK_START", "STAG_PEAK_END", "STAG_PEAK_LENGTH"] + [e+"_MOTIFS" for e in motifLabelList]
if(headerRep == 2):
  header1 = ["LOOP1_"+e for e in header1] + ["LOOP2_"+e for e in header1]
  header2 = ["\t".join(["LOOP1_"+k for k in e.split("\t")]) for e in header2] + ["\t".join(["LOOP2_"+k for k in e.split("\t")]) for e in header2]
  header3 = ["LOOP1_"+e for e in header3] + ["LOOP2_"+e for e in header3]

# Iterating in main file
stagFile = open(stagFileName, "rU")
outputFile1 = open(outputFileName1, "w")
outputFile1.write("\t".join(header1)+"\n")
outputFile2 = open(outputFileName2, "w")
outputFile2.write("\t".join(header2)+"\n")
outputFile3 = open(outputFileName3, "w")
outputFile3.write("\t".join(header3)+"\n")
for line in stagFile:

  # Initialization
  vectorTable1 = []
  vectorTable2 = []
  vectorTable3 = []
  ll = line.strip().split("\t")
  if(fileIsLoop == "Y"): coordinateList = [ll[:3], ll[3:]]
  else: coordinateList = [ll[:3]]
  if(coordinateList[0][0] not in chrList): continue

  # Iterating on coordinates
  for coordVec in coordinateList:

    ###################################################################################################
    # Table Genome
    ###################################################################################################

    # Basic Coordinates
    chrom = coordVec[0]; vectorTable1.append(chrom); vectorTable2.append(chrom); vectorTable3.append(chrom)
    start = int(coordVec[1]); vectorTable1.append(str(start)); vectorTable2.append(str(start)); vectorTable3.append(str(start))
    end = int(coordVec[2]); vectorTable1.append(str(end)); vectorTable2.append(str(end)); vectorTable3.append(str(end))
    regionLength = end - start; vectorTable1.append(str(regionLength)); vectorTable2.append(str(regionLength)); vectorTable3.append(str(regionLength))

    # %CG content
    sequence = str(genomeFile.fetch(chrom, start, end)).upper()
    cgFreq = 0
    for character in sequence:
      if(character == "C" or character == "G"): cgFreq += 1
    vectorTable1.append(str(round(float(cgFreq)/regionLength,2)))

    # CTCF status
    higherScore = -999.
    higherStrand = "NA"
    for k in ctcfIndexList:
      motifFile = motifFileList[k]
      motifFetch = motifFile.fetch(chrom, start, end)
      for read in motifFetch:
        rr = read.qname.split(":")
        motifScore = float(rr[1])
        if(motifScore > higherScore):
          higherScore = motifScore
          higherStrand = "+"
          if(read.is_reverse): higherStrand = "-"
    vectorTable1.append(higherStrand)
    
    # Genomic region
    regionVec = []
    regionFetch = regionsFile.fetch(chrom, start, end)
    for read in regionFetch:
      qName = "NA"
      if(read.qname and read.qname != "."): qName = read.qname
      try: geneSymbol = aliasDict[read.qname.split(":")[1]]
      except Exception: geneSymbol = "NA"
      if(qName == "INTERGENIC:."): qName = "INTERGENIC:NA"
      startx = "NA"
      if(read.pos and read.pos != "."): startx = str(read.pos)
      endx = "NA"
      if(read.aend and read.aend != "."): endx = str(read.aend)
      strand = "+"
      if(read.is_reverse): strand = "-"
      if(qName == "INTERGENIC:NA"): strand = "NA"
      regionVec.append(":".join([qName, geneSymbol, startx, endx, strand]))
    if(regionVec): vectorTable1.append(",".join(regionVec))
    else: vectorTable1.append("NA")

    # ChromHmm Location
    """
    chromHmmStates = []
    chromHmmFetch = chrommHmmFile.fetch(chrom, start, end)
    for read in chromHmmFetch:
      qName = "NA"
      if(read.qname and read.qname != "."): qName = read.qname
      startx = "NA"
      if(read.pos and read.pos != "."): startx = str(read.pos)
      endx = "NA"
      if(read.aend and read.aend != "."): endx = str(read.aend)
      chromHmmStates.append(":".join([qName, startx, endx]))
    if(chromHmmStates): vectorTable1.append(",".join(chromHmmStates))
    else:
    """
    vectorTable1.append("NA")

    # Closest genes
    geneList = []
    geneFetch = regionsFile.fetch(chrom, start - geneExt, end + geneExt)
    for read in geneFetch:
      rr = read.qname.split(":")
      if(rr[0] != "GENE"): continue
      qName = "NA"
      if(rr[-1] and rr[-1] != "."): qName = rr[-1]
      try: geneSymbol = aliasDict[rr[1]]
      except Exception: geneSymbol = "NA"
      startx = "NA"
      if(read.pos and read.pos != "."): startx = str(read.pos)
      endx = "NA"
      if(read.aend and read.aend != "."): endx = str(read.aend)
      strand = "+"
      if(read.is_reverse): strand = "-"
      geneList.append(":".join([qName, geneSymbol, startx, endx, strand]))
    if(geneList): vectorTable1.append(",".join(geneList))
    else: vectorTable1.append("NA")

    # Distance to gene TSS
    gdistVec = []
    center = (start + end) / 2
    for gene in geneList:
      gg = gene.split(":")
      if(gg[4] == "+"): tss = int(gg[2])
      elif(gg[4] == "-"): tss = int(gg[3])
      else:
        gdistVec.append("NA")
        continue
      gdistVec.append(str(tss-center))
    if(gdistVec): vectorTable1.append(",".join(gdistVec))
    else: vectorTable1.append("NA")

    # Stag1/2 expression of genes ([STAG1->[log2FC, minusAux, plusAux], STAG2->[log2FC, minusAux, plusAux]])
    expressionMat = [[] for e in range(0,len(expressionLabelList))]
    for gene in geneList:
      ensg = gene.split(":")[0]
      for k in range(0,len(exprDictList)):
        try:
          exprVec = exprDictList[k][ensg]
          exprVec0 = "NA"
          if(exprVec[0] and exprVec[0] != "."): exprVec0 = exprVec[0]
          exprVec1 = "NA"
          if(exprVec[1] and exprVec[1] != "."): exprVec1 = exprVec[1]
          exprVec2 = "NA"
          if(exprVec[2] and exprVec[2] != "."): exprVec2 = exprVec[2]
        except Exception:
          expressionMat[k].append(":".join(["NA","NA","NA"]))
          continue
        expressionMat[k].append(":".join([exprVec0, exprVec1, exprVec2]))
    for expressionVec in expressionMat:
      if(expressionVec): vectorTable1.append(",".join(expressionVec))
      else: vectorTable1.append("NA")

    # Closest enhancers
    enhancerVec = []
    enhancerFetch = enhancersFile.fetch(chrom, start - geneExt, end + geneExt)
    for read in enhancerFetch:
      nameVec = read.qname.split(":")
      nameVec0 = "NA"
      if(nameVec[0] and nameVec[0] != "."): nameVec0 = nameVec[0]
      startx = "NA"
      if(read.pos and read.pos != "."): startx = str(read.pos)
      endx = "NA"
      if(read.aend and read.aend != "."): endx = str(read.aend)
      nameVec2 = "NA"
      if(nameVec[2] and nameVec[2] != "."): nameVec2 = nameVec[2]
      nameVec1 = "NA"
      if(nameVec[1] and nameVec[1] != "."):
        try: nameVec1 = ensemblDict[nameVec[1]]
        except Exception: nameVec1 = nameVec[1]
      try: geneSymbol = aliasDict[nameVec1]
      except Exception: geneSymbol = "NA"
      nameVec3 = "NA"
      if(nameVec[3] and nameVec[3] != "."): nameVec3 = nameVec[3]
      enhancerVec.append(":".join([nameVec0, startx, endx, nameVec2, nameVec1, geneSymbol, nameVec3]))
    if(enhancerVec): vectorTable1.append(",".join(enhancerVec))
    else: vectorTable1.append("NA")

    # Distance to enhancer center
    enhgenDistVec = []
    for enh in enhancerVec:
      ee = enh.split(":")
      try: enhancerCenter = (int(ee[1]) + int(ee[2])) / 2
      except Exception:
        enhgenDistVec.append("NA")
        continue
      distanceToEnhancerCenter = str(enhancerCenter - center)
      enhgenDistVec.append(distanceToEnhancerCenter)
    if(enhgenDistVec): vectorTable1.append(",".join(enhgenDistVec))
    else: vectorTable1.append("NA")

    # Stag1/2 expression on enhancer targets ([STAG1->[log2FC, minusAux, plusAux], STAG2->[log2FC, minusAux, plusAux]])
    expressionTgtMat = [[] for e in range(0,len(expressionLabelList))]
    for enh in enhancerVec:
      targetGene = enh.split(":")[4]
      for k in range(0,len(exprDictList)):
        try:
          exprVec = exprDictList[k][ensg]
          exprVec0 = "NA"
          if(exprVec[0] and exprVec[0] != "."): exprVec0 = exprVec[0]
          exprVec1 = "NA"
          if(exprVec[1] and exprVec[1] != "."): exprVec1 = exprVec[1]
          exprVec2 = "NA"
          if(exprVec[2] and exprVec[2] != "."): exprVec2 = exprVec[2]
        except Exception:
          expressionTgtMat[k].append(":".join(["NA","NA","NA"]))
          continue
        expressionTgtMat[k].append(":".join([exprVec0, exprVec1, exprVec2]))
    for expressionVec in expressionTgtMat: 
      if(expressionVec): vectorTable1.append(",".join(expressionVec))
      else: vectorTable1.append("NA")

    # GO Terms
    goTermVec = []
    for gene in geneList:
      ensg = gene.split(":")[0]
      try: goTermVec.append(goDict[ensg])
      except Exception: goTermVec.append("NA")
    if(goTermVec): vectorTable1.append(",".join(goTermVec))
    else: vectorTable1.append("NA")

    ###################################################################################################
    # Table Signal
    ###################################################################################################

    # Signals (DNase-seq, histone modifications, TFs) [total signal (RPM), average signal (RPKM)]
    for i in range(0, len(signalFileList)):
      signalFile = signalFileList[i]
      currRpm = rpmList[i]
      try: totalSignal = fetchTotalReadsBam(signalFile, [chrom, start, end])
      except Exception:
        vectorTable2.append("NA,NA")
        continue
      try: vectorTable2.append(str(round(float(totalSignal)/currRpm,2)))
      except Exception: vectorTable2.append("NA")
      try: vectorTable2.append(str(round((float(totalSignal)/currRpm)/regionLength,2)))
      except Exception: vectorTable2.append("NA")

    # Control Signals (DNase-seq, histone modifications, TFs) [total signal (RPM), average signal (RPKM)]
    for i in range(0, len(controlFileList)):
      signalFile = controlFileList[i]
      currRpm = rpmControlList[i]
      try: totalSignal = fetchTotalReadsBam(signalFile, [chrom, start, end])
      except Exception:
        vectorTable2.append("NA,NA")
        continue
      try: vectorTable2.append(str(round(float(totalSignal)/currRpm,2)))
      except Exception: vectorTable2.append("NA")
      try: vectorTable2.append(str(round((float(totalSignal)/currRpm)/regionLength,2)))
      except Exception: vectorTable2.append("NA")

    ###################################################################################################
    # Table Motives
    ###################################################################################################

    # Motives associated to TFs inside region
    for motifFile in motifFileList:
      tfMotifList = []
      try: motifFetch = motifFile.fetch(chrom, start, end)
      except Exception:
        tfMotifList.append("NA")
        continue
      for read in motifFetch:
        rr = read.qname.split(":")
        motifName = "NA"
        if(rr[0] and rr[0] != "."): motifName = rr[0]
        startx = "NA"
        if(read.pos and read.pos != "."): startx = str(read.pos)
        endx = "NA"
        if(read.aend and read.aend != "."): endx = str(read.aend)
        motifScore = "NA"
        if(rr[1] and rr[1] != "."): motifScore = rr[1]
        strand = "+"
        if(read.is_reverse): strand = "-"
        tfMotifList.append(":".join([motifName, startx, endx, motifScore, strand]))
      if(tfMotifList): vectorTable3.append(",".join(tfMotifList))
      else: vectorTable3.append("NA")

  # Writing newline
  outputFile1.write("\t".join(vectorTable1)+"\n")
  outputFile2.write("\t".join(vectorTable2)+"\n")
  outputFile3.write("\t".join(vectorTable3)+"\n")

###################################################################################################
# Table Peak Signal
###################################################################################################

# Signals (DNase-seq, histone modifications, TFs) [total signal (RPM), average signal (RPKM)]
vectorTable4 = []
vectorTable5 = []
for i in range(0, len(signalFileList)):
  vec4 = []
  vec5 = []
  peakFileName = peakFileNameList[i]
  signalFile = signalFileList[i]
  currRpm = rpmList[i]
  totalSignal = 0
  peakRegionLength = 0
  peakFile = open(peakFileName, "rU")
  for line in peakFile:
    ll = line.strip().split("\t")
    region = [ll[0], int(ll[1]), int(ll[2])]
    peakRegionLength = (int(ll[2])-int(ll[1]))
    try: totalSignal = fetchTotalReadsBam(signalFile, region)
    except Exception:
      vec4.append("NA")
      vec5.append("NA")
      continue
    try: vec4.append(str(round(float(totalSignal)/currRpm,2)))
    except Exception: vec4.append("NA")
    try: vec5.append(str(round((float(totalSignal)/currRpm)/peakRegionLength,2)))
    except Exception: vec5.append("NA")
  peakFile.close()
  vectorTable4.append(vec4)
  vectorTable5.append(vec5)

# Writing table peak signal
maxV = max([len(e) for e in vectorTable4])
outputFile4 = open(outputFileName4, "w")
header4 = ["TOTAL_PEAK_"+e+"_NORMALIZED_COUNT(RPM)\t"+"AVERAGE_PEAK_"+e+"_NORMALIZED_COUNT(RPKM)" for e in signalLabelList]
outputFile4.write("\t".join(header4)+"\n")
for i in range(0,maxV):
  vec = []
  for j in range(0,len(vectorTable4)):
    try: vec.append(vectorTable4[j][i])
    except Exception: vec.append("NA")
    try: vec.append(vectorTable5[j][i])
    except Exception: vec.append("NA")
  outputFile4.write("\t".join(vec)+"\n")

stagFile.close()
outputFile1.close()
outputFile2.close()
outputFile3.close()
outputFile4.close()
genomeFile.close()
regionsFile.close()
#chrommHmmFile.close()
enhancersFile.close()
[e.close() for e in signalFileList]
[e.close() for e in controlFileList]
[e.close() for e in motifFileList]

# Removing all files
command = "rm -rf "+tempLocation
os.system(command)


