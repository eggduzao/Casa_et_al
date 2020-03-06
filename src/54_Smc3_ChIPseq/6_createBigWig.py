
#Import
import os
import sys
from pysam import Samfile
from pysam import __version__ as ps_version

###################################################################################################
# Input
###################################################################################################

# Input
downstreamExt = int(sys.argv[1])
upstreamExt = int(sys.argv[2])
genomeSizesFileName = sys.argv[3]
bamFileName = sys.argv[4]
bwFileName = sys.argv[5]

###################################################################################################
# Functions
###################################################################################################

class PileupRegion:

  def __init__(self, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift):
    self.start = start
    self.end = end
    self.length = end - start
    self.downstream_ext = downstream_ext
    self.upstream_ext = upstream_ext
    self.forward_shift = forward_shift
    self.reverse_shift = reverse_shift
    self.vector = [0.0] * self.length

  def __call__(self, alignment):
    try:
      if(not alignment.is_reverse):
        for i in range(max(alignment.pos + self.forward_shift - self.upstream_ext, self.start),  min(alignment.pos + self.forward_shift + self.downstream_ext, self.end - 1)): self.vector[i - self.start] += 1.0
      else:
        for i in range(max(alignment.aend + self.reverse_shift - self.downstream_ext, self.start), min(alignment.aend + self.reverse_shift + self.upstream_ext, self.end - 1)): self.vector[i - self.start] += 1.0
    except Exception: pass

def create_wig(pysamFile, downstream_ext, upstream_ext, forward_shift, reverse_shift, genome_sizes_file_name, output_file_name, window_size = 1, wig_type = "fixedStep"):

  # Fixed parameters
  MAX_GENOME_REGION_SIZE = 1000000

  # Reading genome sizes as a reference
  genome_sizes_dict = dict()
  genome_sizes_file = open(genome_sizes_file_name, "rU")
  for line in genome_sizes_file:
    ll = line.strip().split("\t")
    genome_sizes_dict[ll[0]] = int(ll[1])
  genome_sizes_file.close()
  chrom_list = sorted(genome_sizes_dict.keys())

  # Creating output file
  output_file = open(output_file_name, "w")

  # Iterating on each chromosome
  for chrom in chrom_list:

    # Iterating on genomic regions for memory purposes
    for i in range(0, genome_sizes_dict[chrom], MAX_GENOME_REGION_SIZE):

      # Wig header
      wig_header = "fixedStep chrom="+chrom+" start="+str(i+1)+" step="+str(window_size)
      output_file.write(wig_header+"\n")

      # Region to fetch the signal
      start = i
      end = min(i+MAX_GENOME_REGION_SIZE, genome_sizes_dict[chrom])

      # Fetch signal
      pileup_region = PileupRegion(start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift)
      if(ps_version == "0.7.5"):
        pysamFile.fetch(reference=chrom, start=start, end=end, callback=pileup_region)
      else:
        iter = pysamFile.fetch(reference=chrom, start=start, end=end)
        for alignment in iter:
          pileup_region.__call__(alignment)

      # Print signal
      for s in pileup_region.vector: output_file.write(str(s)+"\n")

  # Closing output file
  output_file.close()

def convert_to_bigwig(wig_file_name, genome_sizes_file_name, bw_file_name, remove_original = False):
  command = " ".join(["wigToBigWig", wig_file_name, genome_sizes_file_name, bw_file_name])
  os.system(command)
  if(remove_original): os.system("rm "+wig_file_name)

def create_bigwig(pysamFile, downstream_ext, upstream_ext, forward_shift, reverse_shift, genome_sizes_file_name, output_file_name, window_size = 1, wig_type = "fixedStep", remove_wig = False):

  # Creating wig file from bam
  temporary_wig_file_name = output_file_name + "_temporary.wig"
  create_wig(pysamFile, downstream_ext, upstream_ext, forward_shift, reverse_shift, genome_sizes_file_name, temporary_wig_file_name, window_size, wig_type)

  # Converting to bigwig
  convert_to_bigwig(temporary_wig_file_name, genome_sizes_file_name, output_file_name, remove_original = remove_wig)

###################################################################################################
# Execution
###################################################################################################

# Initialization
outLoc = "/".join(bwFileName.split("/")[:-1]) + "/"
command = "mkdir -p "+outLoc
os.system(command)

# Parameters
forwardShift = 0
reverseShift = 0

# Creating bigwig
bamFile = Samfile(bamFileName, "rb")
create_bigwig(bamFile, downstreamExt, upstreamExt, forwardShift, reverseShift, genomeSizesFileName, bwFileName, window_size = 1, wig_type = "fixedStep", remove_wig = True)
bamFile.close()


