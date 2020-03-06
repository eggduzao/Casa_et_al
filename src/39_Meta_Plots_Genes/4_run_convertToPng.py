
# Import
import os
import sys
from glob import glob

# Input
counter = 0
pdfFileNameList = glob("/home/egg/Desktop/2_Meta_Gene_Plots_PNG/*/*.pdf")

# Iterating on pdfFileNameList
for pdfFileName in pdfFileNameList:

  if(os.path.isfile(pdfFileName)): counter += 1

  pngFileName = pdfFileName[:-3] + "png"
  command = "pdftoppm "+pdfFileName+" "+pngFileName+" -png"
  os.system(command)

print counter

# Input
counter = 0
pngFileNameList = glob("/home/egg/Desktop/2_Meta_Gene_Plots_PNG/*/*.png")

for pngFileName in pngFileNameList:

  if(os.path.isfile(pngFileName)): counter += 1

  newPngFileName = pngFileName.split("-1")[0]
  command = "mv "+pngFileName+" "+newPngFileName
  os.system(command)

print counter


