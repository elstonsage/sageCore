#============================================================================
# File:      phenotype1.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/15/4  Created
#            10/25/4 Modified phenotype0.py to create this more efficient
#                    version.
#                                                                          
# Notes:     For pooled version for haplotype freq. estimating program.
#            Has much more efficient algorithm for calculating haplotype
#            combinations that are consistent with pool phenotypes.
#
#============================================================================

import random
from copy import copy, deepcopy
from haplotype import *
from count import count
from combinations import combinations
from permute import permute

# - Represents a phenotype including it's constituent haplotype sets and their
#   weights as used by the EM algorithm.  Constructor takes a sequence of 
#   genotype sequences.  'None' represents missing 
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
    phased_sets = []
    for genotype in pheno_seq:
      phases = []
      
      new_genotype = phenotype.convert_genotype(genotype)
      permute(new_genotype, 0, phases)
      phased_sets.append(copy(phases))
      
    self.generate(phased_sets, hap_seq_sets)
    
  # - Convert genotype from a list of allele counts to a list
  #   of allele indices.
  def convert_genotype(genotype):
    new_genotype = []
    for a in range(len(genotype)):
      for n in range(genotype[a]):
        new_genotype.append(a)
        
    return  new_genotype
  
  convert_genotype = staticmethod(convert_genotype)    
      
  # - Generate haplotypes recursively by breaking first set of genotypes
  #   encountered into a set with one genotype and a set with the remaining
  #   genotypes.  Repeat the process until there is one genotype at
  #   every locus, ie, the haplotypes are unambiguous.
  #
  def generate(self, phased_sets, hap_seq_sets):
    for l in range(len(phased_sets)):
      while len(phased_sets[l]) > 1:
        new_sets = deepcopy(phased_sets)
        new_sets[l] = deepcopy([phased_sets[l][-1]])
        
        phased_sets[l] = phased_sets[l][:-1]
        
        self.generate(new_sets, hap_seq_sets)

    new_set = self.make_set(phased_sets)
    new_set.sort()
    if new_set not in hap_seq_sets:
      hap_seq_sets.append(new_set)        
    
  # - From a list of lists of phased genotypes where each sublist is of size 1,
  #   construct a list of haplotypes.
  #
  def make_set(self, phased_sets):
    hap_set = []  
    for h in range(len(phased_sets[0][0])):
      hap_set.append([])
      for l in range(len(phased_sets)):
        hap_set[h].append(phased_sets[l][0][h])
        
    for s in range(len(hap_set)):
      hap_set[s] = tuple(hap_set[s])
        
    return  hap_set
      
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

  pheno_seq = [(3,1), (3,1)]
  pheno = phenotype(pheno_seq, haps, loci)
  
  """
  haps[(1, 1)].new_count = 20
  haps[(0, 0)].new_count = 10    
  haps[(1, 0)].new_count = 5
  """
  
  haps[(0, 1)].new_count = 10
  haps[(0, 0)].new_count = 10
  haps[(1, 1)].new_count = 10
  haps[(1, 0)].new_count = 10 
  
  haps.output()
  pheno.output()
  pheno.calc_weights(sample_size)
  pheno.output()
  