#============================================================================
# File:      phenotype2a.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/15/4  Created
#            10/25/4 Modified phenotype0.py to create this more efficient
#                    version.
#            10/26/4 Modified to make factorial calculation more efficient.
#
#            11/19/4  Modified phenotype1b to incorporate new method of en-
#                     umerating haplotypes sets which was proposed by
#                     Katrina (see generate2.py).
#            11/22/4  Now uses estimate2a.py
#
#                                                                          
# Notes:     For pooled version for haplotype freq. estimating program.
#            Has much faster algorithm for determining haplotype combinations
#            which are consistent with pool phenotypes.
#
#============================================================================

import random
from copy import copy, deepcopy
from haplotype import *
from count import count
from combinations import combinations
from permute import permute
import generate2a

factorials = { 0:1, 1:1, 2:2 }

# - Represents a phenotype including it's constituent haplotype sets and their
#   weights as used by the EM algorithm.  Constructor takes a sequence of 
#   tuples representing allele counts at each locus.  'None' represents missing 
#   information.  Master dictionary of haplotypes is built as phenotypes are 
#   constructed. 
#
class phenotype:
  phenotype_count = 0

  def __init__(self, pheno_seq, haps, loci):
    phenotype.phenotype_count += 1
    print 'creating phenotype %d ...' % phenotype.phenotype_count
    self.count = 1
    self.loci = loci
    self.g_size = self.loci[0].g_size
    self.haps = haps
    self.name = phenotype.create_name(pheno_seq, loci)
    self.hap_sets = {}      # (hap_seq1, hap_seq2 ... hap_seqg_size) : weight
    self.hap_set_counts = {}   # (hap_seq1, hap_seq2 ... hap_seqg_size) : dictionary of haplotype counts
    
    
    hap_seq_sets = []    # keys to hap_sets
    self.generate_hap_seq_sets(pheno_seq, hap_seq_sets)
    
    # - Add any new haplotypes to the master haplotype dictionary.
    #
    for set in hap_seq_sets:
      for hap_seq in set:
        if not hap_seq in self.haps.haps:
          self.haps[hap_seq] = haplotype(hap_seq, self.loci)
        
      self.hap_sets[tuple(set)] = None        
        
    self.ambiguous = len(self.hap_sets) > 1
    self.init_weights()
    
    for set in self.hap_sets.keys():
      self.hap_set_counts[set] = phenotype.copy_counts(set)    
    
  def create_name(pheno_seq, loci):
    name = ""
    for l in range(len(pheno_seq)):
      name += loci[l].g_create_name(pheno_seq[l])
      name += '  '
      
    return  name
    
  create_name = staticmethod(create_name)
    
  def __repr__(self):
    return  '%d\t%s' % (self.count, self.hap_sets.__repr__())

  # - Generate all possible haplotype sets consistent with the phenotype
  #   sequence.
  #
  def generate_hap_seq_sets(self, pheno_seq, hap_seq_sets):
    generate2a.generate(pheno_seq, hap_seq_sets)
      
  def output(self):
    print self.name
    print '\ncount   %f' % self.count
    print 'haplotype sets'
    for set in self.hap_sets:
      print
      for h in set:
        print h,
        print '  %s' % self.haps[h].name
        
      print 'weight %f' % self.hap_sets[set] 
  
  
  # - Assign random weights with which to begin EM algorithm.
  #  
  def init_weights(self):
    if self.ambiguous:
      set_count = len(self.hap_sets)
      weights = []
      total_weight = 0.0
      for i in range(set_count):
        weight = 1 - random.uniform(0, 1)
        weights.append(weight)
        total_weight += weight
        
      keys = self.hap_sets.keys()
      for k in range(set_count):
        self.hap_sets[keys[k]] = weights[k] / total_weight
        
    else:   # Only one set.
      keys = self.hap_sets.keys()
      self.hap_sets[keys[0]] = 1.0
      
      
  def incr_count(self):
    self.count += 1      
      
  # - Determine weights as a function of estimated haplotype
  #   frequencies.
  #
  def calc_weights(self, sample_size):
    """ Ito et. al.  AJHG 72:384-398, 2003
        Equation (1)  """
    
    if self.ambiguous:
      for set in self.hap_sets:
        self.hap_sets[set] = self.set_prob(set, self.g_size, sample_size)
        
      # - Normalize the weights.
      #
      prob_total = reduce(lambda p1, p2 : p1 + p2, self.hap_sets.values())
      for set in self.hap_sets:
        self.hap_sets[set] /= prob_total
        
  # - Determine probability of an individual haplotype combination.
  #
  def set_prob(self, set, g_size, sample_size):
    """ Ito et. al.  AJHG 72:384-398, 2003  
        Step 3 of EM Algorithm section w. '/' omitted!  """
    
    counts = self.hap_set_counts[set]
    distinct_haps = len(self.hap_set_counts[set])
    
    product = 1
    for hap in counts:
      product *= self.haps[hap].new_freq(sample_size) ** counts[hap] \
                 / phenotype.factorial(counts[hap])
      
    return  phenotype.factorial(g_size) * product
    
  def factorial(number):
    if factorials.has_key(number):
      return  factorials[number]
    else:
      factorials[number] = phenotype.factorial(number - 1) * number
      return  phenotype.factorial(number)
    
  factorial = staticmethod(factorial)
  
  # - Return a dictionary whose keys are haplotypes in the 'set' and whose
  #   values are the number of copies of the haplotype in the set.
  #
  def copy_counts(set):
    counts = {}
    for h in set:
      if counts.has_key(h):
        counts[h] += 1
      else:
        counts[h] = 1
        
    return  counts
  
  copy_counts = staticmethod(copy_counts)


if __name__ == '__main__':
  """
  sample_size = 10
  pool_size = 2
  
  loci = []
  
  loc = locus(pool_size)
  loc.add_allele('A')
  loc.add_allele('a')
  
  loci.append(loc)
  loci.append(loc)
  
  haps = haplotype_set(sample_size)

  pheno_seq = [(3,1), (3,1)]
  pheno = phenotype(pheno_seq, haps, loci)
  
  haps[(1, 1)].new_count = 20
  haps[(0, 0)].new_count = 10    
  haps[(1, 0)].new_count = 5
  
  haps[(0, 1)].new_count = 10
  haps[(0, 0)].new_count = 10
  haps[(1, 1)].new_count = 10
  haps[(1, 0)].new_count = 10 
  
  haps.output()
  pheno.output()
  pheno.calc_weights(sample_size)
  pheno.output()
  """
  
  for i in range(7):
    print '%d! = %d' % (i, phenotype.factorial(i))