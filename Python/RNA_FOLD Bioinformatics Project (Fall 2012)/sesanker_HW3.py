# This script identifies all noncoding RNAs by using a sliding Z-score
# Include sequences files in same folder as script. 
# Sequence files should be in FASTA format with each sequence seperated by a '>' character
# Each sequence is referred to as a candidate 'Region'
# to find all ncRNAs use the function: 
# find_all_ncRNA(["file 1",...,"file n"],[(number sequences in file 1) + 1,...,(number sequences in file n) + 1])
# see bottom of script for the 5 files examined for this assignment
# Outputs a file for each input file. In each file is a list of ncRNAs for each Region in the file

# Atempted to complete in as few lines of code as possible. Also, no loops and no global variables are used. Only functional constructs

# Colbert Sesanker 7/12

from __future__ import print_function
from subprocess import Popen, PIPE, STDOUT
import re, random 
import ushuffle as u

def calculate_MFE(seq):   
  scpExec = Popen("RNAfold", stdin=PIPE, stdout=PIPE, stderr = PIPE, close_fds = True)   
  stdout, stderr = scpExec.communicate(seq) 
  mean_folding_energy = re.search('-?\d*\.\d+', stdout ).group(0)  
  return float(mean_folding_energy) 

def z_score(seq):
  shuffles = 500
  def add(x,y): return x + y  
  shuffled_seqs = [''.join(di.dinuclShuffle(seq)) for i in range(shuffles)]  
  mfes = map(calculate_MFE, shuffled_seqs)
  mean_mfe = reduce(add, mfes, 0) / shuffles
  sd_mfe = (reduce(add, map(lambda x: (x - mean_mfe)**2, mfes), 0) / shuffles)**(.5)
  return (calculate_MFE(seq) - mean_mfe) / sd_mfe
  
def find_ncRNA(seq, num):
  windows = [seq[j:i+j] for i in [80, 120, 160, 200] for j in range(0,len(seq)-i, 5)]  
  window_scores = map(z_score, windows)  
  ncRNA = [[window_scores[i], windows[i]] for i in range(len(windows)) if window_scores[i] < -3.5]
  return ["Region " + str(num) + ":"] + ncRNA

def find_ncRNAs(seq_file, end):
   seqs = map(lambda x: re.sub('\n','',x), re.findall('\n([\S|\s]*?)\n\>', open(seq_file, 'r').read() + "\n>"))
   return map(find_ncRNA, seqs, range(1, end))
   
def parse_ncRNAs(ncRNAs):  
   r = re.sub('\[\'[R]', '\n\n\nR', re.sub('\:\'\,',":\n", re.sub('\]\,',"\n", str(ncRNAs))))
   return re.sub('[\'|\[|\]]', '', r)

def find_all_ncRNA(seq_files, ranges):     
  parsed_ncRNAs = map(parse_ncRNAs, map(find_ncRNAs, seq_files, ranges))
  write = [open(species + "_ncRNAs.txt", 'w') for species in seq_files]
  [print(parsed_ncRNAs[i], file = write[i]) for i in range(len(seq_files))]

find_all_ncRNA(["positive_controls","negative_controls"],[4,2])
find_all_ncRNA(["misc_rna_set11_Sesanker_Colbert","Gatti-Syntenic-seqs.txt","Neoformans-JEC21-Syntenic-seqs.txt"],[21,8,8])

