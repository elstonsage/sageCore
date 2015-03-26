#============================================================================
# File:      phenotype.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   12/4/3  Created
#            12/18/3 Modified to allow an arbitrary number of alleles at a locus.
#                                                                          
# Notes:   
#
#============================================================================

import random
from copy import copy
from haplotype import *

# - Represents a phenotype including it's constituent haplotype pairs and their
#   weights as used by the EM algorithm.  Constructor takes a sequence of 
#   genotype indices.  'None' represents missing information.  External 
#   dictionary of haplotypes is built as phenotypes are constructed. 
#
class phenotype:

  # - Generate all possible haplotype pairs corresponding to the phenotype
  #   sequence.
  #
  def generate_hap_seq_pairs(self, pheno_seq, hap_seq_pairs, loci):
  
    # - Split sequence into multiple sequences at each locus with missing info.
    #
    m = self.find_missing(pheno_seq)
    if m != -1:     # Sequence contains missing information.
      for g in range(len(loci[m].genos)):
        new_pheno_seq = copy(pheno_seq)
        new_pheno_seq[m] = g
        self.generate_hap_seq_pairs(new_pheno_seq, hap_seq_pairs, loci)
    
    else:           # No missing information.
      
      # - Find the number and location of the heterozygotes.  Use a binary
      #   counting method to generate haplotypes w. each possible combination
      #   of phases.
      #
      h_loci = self.find_heterozygotes(pheno_seq, loci)
      h_count = len(h_loci)
      
      if h_count > 0:   
      
        # - Generate sequences for all permutations of phase.
        #
        for i in range(2 ** (h_count - 1)):
          pheno_seq0 = copy(pheno_seq)
          pheno_seq1 = copy(pheno_seq)
          hap_seq0 = self.generate_hap_seq(pheno_seq0, h_loci, i, loci)
          hap_seq1 = self.generate_hap_seq(pheno_seq1, h_loci, 2 ** h_count - 1 - i, loci)
          
          hap_seq_pairs.append([hap_seq0, hap_seq1])
          
      else:
        hap_seq = self.generate_hap_seq(copy(pheno_seq), h_loci, None, loci)
        hap_seq_pairs.append([hap_seq, hap_seq])
           
        
  # - Generate haplotype sequence corresponding to given phase information.
  #
  def generate_hap_seq(self, pheno_seq, h_loci, phase_info, loci):
    hap_seq = pheno_seq
    if phase_info == None:    # Homozygous at every locus.
      phases = []
    else:
      phases = int_to_seq(phase_info, len(h_loci))
      
    phase_index = 0
    for l in range(len(pheno_seq)):
      if loci[l].genos[pheno_seq[l]].homozygous: 
        hap_seq[l] = loci[l].genos[pheno_seq[l]][0]
      else:
        hap_seq[l] = loci[l].genos[pheno_seq[l]][phases[phase_index]]
        phase_index += 1                  

    return  hap_seq

  # - Return location of first missing value.  If none are found,
  #   return -1
  #
  def find_missing(self, pheno_seq):
    m = -1
    for l in range(len(pheno_seq)):
      if pheno_seq[l] == None:
        m = l
        break
        
    return  m
  
  # - Find locations of all heterozygotes.
  #
  def find_heterozygotes(self, pheno_seq, loci):
    locs = []
    for l in range(len(pheno_seq)):
      if loci[l].genos[pheno_seq[l]].heterozygous:
        locs.append(l)

    return  locs

  def __init__(self, pheno_seq, haps, loci):
    self.count = 1
    self.hap_pairs = {}
    
    hap_seq_pairs = []
    self.generate_hap_seq_pairs(pheno_seq, hap_seq_pairs, loci)
    
    # - Add any new haplotypes to set of haplotypes.
    #
    for pair in hap_seq_pairs:
      first = tuple(pair[0])
      if not haps.has_key(first):
        haps[first] = haplotype(first, loci)
      second = tuple(pair[1])
      if not haps.has_key(second):
        haps[second] = haplotype(second, loci)
        
      self.hap_pairs[(first, second)] = None        
        
    self.ambiguous = len(self.hap_pairs) > 1
    self.init_weights()
    
  def __repr__(self):
    return  '%d\t%s' % (self.count, self.hap_pairs.__repr__())

  def output(self, haps):
    print 'count   %f' % self.count
    print 'haplotype pairs'
    for key in self.hap_pairs:
      print '%s  %s   %f' % (haps[key[0]].name, haps[key[1]].name, self.hap_pairs[key]) 
  
  
  # - Assign random weights with which to begin EM algorithm.
  #  
  def init_weights(self):
    if self.ambiguous:
      pair_count = len(self.hap_pairs)
      weights = []
      total_weight = 0.0
      for i in range(pair_count):
        weight = 1 - random.uniform(0, 1)
        weights.append(weight)
        total_weight += weight
        
      keys = self.hap_pairs.keys()
      for k in range(pair_count):
        self.hap_pairs[keys[k]] = weights[k] / total_weight
        
    else:   # Only one pair.
      keys = self.hap_pairs.keys()
      self.hap_pairs[keys[0]] = 1.0
      
      
  def incr_count(self):
    self.count += 1      
      
  # - Determine weights as a function of estimated haplotype
  #   frequencies.
  #
  def calc_weights(self, haps, sample_size):
    if self.ambiguous:
      total_freq = 0.0
      for key in self.hap_pairs.keys():
        pair_term = haps[key[0]].new_freq(sample_size) * \
                    haps[key[1]].new_freq(sample_size)
        if key[0] != key[1]:
          pair_term *= 2
        total_freq += pair_term
        
      for key in self.hap_pairs.keys():
        pair_term = haps[key[0]].new_freq(sample_size) * \
                    haps[key[1]].new_freq(sample_size)
        if key[0] != key[1]:
          pair_term *= 2      
        self.hap_pairs[key] = pair_term / total_freq


def pheno_seq_to_string(pheno_seq, loci):
  str = ""
  for l in range(len(pheno_seq)):
    if pheno_seq[l] == None:
      str += '?/? '
    else:
      str += '%s ' % loci[l].g_name(pheno_seq[l])
  
  return  str


if __name__ == '__main__':
  
  loci = []
  
  loc = locus()
  loc.add_allele('A')
  loc.add_allele('a')
  loc.build_genotypes()
  
  loc2 = locus()
  loc2.add_allele('a')
  loc2.add_allele('A')
  
  loci.append(loc)
  loci.append(loc)
  loci.append(loc)
  loci.append(loc)
  
  haps = haplotype_set(10)
  pheno_seq = [None, None, None]
  pheno = phenotype(pheno_seq, haps, loci)
  print haps
  print pheno
  