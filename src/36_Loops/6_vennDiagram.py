
###################################################################################################
# Input
###################################################################################################

# Import
import os
import sys
sys.path = ["/home/egusmao/.local/lib/python2.7/site-packages"] + sys.path
from pysam import Samfile

# Input
resolution = int(sys.argv[1])
loopFileName = sys.argv[2]
stagRegionList = sys.argv[3].split(",")
tempLocation = sys.argv[4]
outputFileName = sys.argv[5]

# Initialization
command = "mkdir -p "+tempLocation
os.system(command)

###################################################################################################
# Functions
###################################################################################################

def check_bam_at_least_one_read(bam_file, region):
  res = False
  for read in bam_file.fetch(region[0], region[1], region[2]):
    res = True
    break
  return res

def convert_loop_file_into_ranges(resolution, loop_file_name, temporary_location):
  
  # Execution
  out_file_name = temporary_location + "out_file_name.bed"
  loop_file = open(loop_file_name, "rU")
  out_file = open(out_file_name, "w")
  halfRes = resolution / 2
  for line in loop_file:
    ll = line.strip().split("\t")
    chrom = ll[0]; p1 = int(ll[1]); p2 = int(ll[2])
    out_file.write("\t".join([str(e) for e in [chrom, max(0, p1 - halfRes), p1 + halfRes, max(0, p2 - halfRes), p2 + halfRes]])+"\n")
  loop_file.close()
  out_file.close()

  # Return objects
  return out_file_name

def calculate_intersection(loop_file_name, stag_region_list):

  # Creating intersection vector
  intersection_vector = [0 for e in stag_region_list] + [0]

  # Opening files
  loop_file = open(loop_file_name, "rU")
  stag_region_file_vector = [Samfile(e,"rb") for e in stag_region_list]

  # Calculating intersections
  for line in loop_file:
    ll = line.strip().split("\t")
    range1 = [ll[0], int(ll[1]), int(ll[2])]
    range2 = [ll[0], int(ll[3]), int(ll[4])]
    flagIntAtLeastOne = False
    for i in range(0, len(stag_region_file_vector)):
      stag_file = stag_region_file_vector[i]
      int1 = check_bam_at_least_one_read(stag_file, range1)
      int2 = check_bam_at_least_one_read(stag_file, range2)
      if(int1 or int2):
        intersection_vector[i] += 1
        flagIntAtLeastOne = True
    if(not flagIntAtLeastOne): intersection_vector[-1] += 1

  # Closing files
  loop_file.close()
  for e in stag_region_file_vector: e.close()

  # Returning objects
  return intersection_vector

def writing_intersection_vector(intersection_vector, output_file_name):

  # Writing intersection to file
  iv = intersection_vector
  header = ["STAG1", "STAG2", "NONE"]
  vec1 = ["1", str(iv[0]), "1000000"]
  vec2 = [str(iv[0]+iv[2]), str(iv[0]+iv[1]+iv[2]), str(1000000 + iv[-1])]
  output_file = open(output_file_name, "w")
  output_file.write("\t".join(header)+"\n")
  output_file.write("\t".join(vec1)+"\n")
  output_file.write("\t".join(vec2)+"\n")
  output_file.close()

def create_venn_digram_table(resolution, loop_file_name, stag_region_list, temporary_location, output_file_name):

  # Calculate how much intersections
  new_loop_file_name = convert_loop_file_into_ranges(resolution, loop_file_name, temporary_location)

  # Calculate intersection
  intersection_vector = calculate_intersection(new_loop_file_name, stag_region_list)

  # Writing intersection_vector
  writing_intersection_vector(intersection_vector, output_file_name)

###################################################################################################
# Execution
###################################################################################################

# Create Matrices
create_venn_digram_table(resolution, loopFileName, stagRegionList, tempLocation, outputFileName)

# Remove temporary folder
command = "rm -rf "+tempLocation
os.system(command)


