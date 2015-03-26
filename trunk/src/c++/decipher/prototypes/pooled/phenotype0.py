#============================================================================
# File:      phenotype0.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/15/4  Created
#                                                                          
# Notes:     For pooled version for haplotype freq. estimating program.
#
#============================================================================

import random
from copy import copy
from haplotype import *
from count import count
from combinations import combinations

# - Represents a phenotype including it's constituent haplotype sets and their
#   weights as used by the EM algorithm.  Constructor takes a sequence of 
#   genotype sequences (as defined in locus.py).  'None' represents missing 
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
    
    hap_seq_sets = []
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
    """ hap_list [[a1,a2,a3],[a1,a2,a3], ... last consistent haplotype]
        set      [h0 count,h1 count ... hn count]
        seq_set  [[a1,a2,a3],[a1,a2,a3] ... g_sizeth haplotype] """
        
    hap_list = self.generate_haps(pheno_seq)
    for set in count(len(hap_list), self.g_size + 1):
      if reduce(lambda c1, c2 : c1 + c2, set) == self.g_size:
        if phenotype.is_consistent(set, hap_list, pheno_seq):
          hap_seq_sets.append(phenotype.set_to_seq_set(set, hap_list))
      
  def is_consistent(set, hap_list, pheno_seq):
    for l in range(len(pheno_seq)):
      for a in range(len(pheno_seq[l])):
        req_count = pheno_seq[l][a]
        act_count = 0
        for h in range(len(set)):
          if hap_list[h][l] == a:
            act_count += set[h]
        if req_count != act_count:
          return  False
          
    return  True
  
  is_consistent = staticmethod(is_consistent)
  
  def set_to_seq_set(set, hap_list):
    seq_set = []
    for h in range(len(set)):
      for c in range(set[h]):
        seq_set.append(tuple(hap_list[h]))
        
    return  seq_set
  
  set_to_seq_set = staticmethod(set_to_seq_set)

  def generate_haps(self, pheno_seq):
    """ pheno_seq example
  
                |<------- 2 alleles --------->|
                
                            |<- g_size = 4 ->|                
                            
                [[1, 2, 1], [0, 0, 1, 1, 1, 1]]
                     ^               
                     |_ 2 copies of allele 1      """
                     
    allele_list = self.generate_allele_list(pheno_seq)
    hap_list = combinations(allele_list)
    return  hap_list
    
  def generate_allele_list(self, pheno_seq):
    """ allele_list example
    
      [(0, 1, 2), (2, 3, 4, 5)] 
      
        indices of alleles present at each locus  """
        
    allele_list = []
    for l in range(len(pheno_seq)):
      alleles = [ph for ph in range(len(pheno_seq[l])) if pheno_seq[l][ph] != 0]
      allele_list.append(tuple(alleles))
      
    return  allele_list

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
    
    counts = phenotype.copy_counts(set)
    distinct_haps = len(counts)
    
    product = 1
    for hap in counts:
      product *= self.haps[hap].new_freq(sample_size) ** counts[hap] \
                 / phenotype.factorial(counts[hap])
      
    return  phenotype.factorial(g_size) * product
    
  def factorial(number):
    return  reduce(lambda n1, n2 : n1 * n2, range(1, number + 1))
    
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
  sample_size = 10
  pool_size = 2
  
  loci = []
  
  loc = locus(pool_size)
  loc.add_allele('A')
  loc.add_allele('a')
  
  loci.append(loc)
  loci.append(loc)
  
  haps = haplotype_set(sample_size)

  pheno_seq = [(1, 3), (3,1)]
  pheno = phenotype(pheno_seq, haps, loci)
  
  haps[(1, 1)].new_count = 20
  haps[(0, 0)].new_count = 10    
  haps[(1, 0)].new_count = 5
  haps[(0, 1)].new_count = 5        
  
  haps.output()
  pheno.output()
  pheno.calc_weights(sample_size)
  pheno.output()
  