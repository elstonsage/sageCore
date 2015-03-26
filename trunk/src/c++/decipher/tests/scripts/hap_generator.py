#============================================================================
# File:      hap_generator.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 5-26-5
#                                                                          
# Notes:     Generate a random list of haplotypes and their frequencies
#            from among all possible haplotypes for a given number of SNP's.
#
#============================================================================

import sys
import random

# - return list of 1's and 2's corresponding to the binary representation of 
#   number.
#
def sequence(number, digit_count):
  seq = []
  
  for d in range(digit_count):
    least_sig = number % 2
    if(least_sig):
      seq.append('2')
    else:
      seq.append('1')
    
    number = number >> 1
    
  seq.reverse()
  
  return seq
  

def create_haps(locus_count, hap_count):
  max_hap_count = 2 ** locus_count
  assert(locus_count <= max_hap_count)

  haps = []
  hap_numbers = []
    
  # - 5 haplotypes account for 90% of the frequency.
  #
  while len(haps) < hap_count:
    if len(haps) < 5:
      frequency = .18
    else:
      frequency = .1 / float(hap_count - 5)
      
    hap_number = random.randint(0, max_hap_count - 1)
    
    if not (hap_number in hap_numbers):
      hap_sequence = sequence(hap_number, locus_count)
      hap_numbers.append(hap_number)
      haps.append([frequency, hap_sequence])

  return  haps  
  
if __name__ == "__main__":
  print create_haps(9, 20)