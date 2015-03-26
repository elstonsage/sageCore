#============================================================================
# File:      haplotype.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/15/4  Created.
#                                                                          
# Notes:     For pooled version for haplotype freq. estimating program.
#
#============================================================================

from locus import *

# - Stores haplotype name and counts.  Constructor takes a sequence of 
#   indices indicating which allele is present at each locus.
#
class haplotype:

  def __init__(self, hap_seq, loci):
    self.loci = loci
    self.name = self.generate_name(hap_seq)
    self.old_count = 0.0
    self.new_count = 0.0
    self.static_count = 0.0   # - Portion of count not changing during 
                              #   algorithm iteration.
    
  def __repr__(self):
    return  '%s\tstatic count   %f\told count   %f\tnew count   %f' % (self.name, self.static_count, 
                                                                       self.old_count, self.new_count)
                                                                       
  # - Generate string representation of a haplotype.
  #
  def generate_name(self, hap_seq):
    str = ""
    for l in range(len(hap_seq)):
      str += self.loci[l].a_name(hap_seq[l])
      
    return  str    
                                
  def iadd(self, value):
    self.new_count += value
    
  def reset(self):
    self.old_count = self.new_count
    self.new_count = self.static_count
    
    
  # - Are new counts sufficiently close to old counts?
  #
  def converged(self, epsilon, sample_size):
    diff = abs(self.new_freq(sample_size) - self.old_freq(sample_size))
    return_value = diff < epsilon
    if not return_value:
      print 'frequency diff = %f' % diff
      
    return  diff < epsilon
    
  def new_freq(self, sample_size):    
    return  float(self.new_count) / (sample_size * self.loci[0].g_size)
    
  def old_freq(self, sample_size):    
    return  float(self.old_count) / (sample_size * self.loci[0].g_size)

 
class haplotype_set:
  def __init__(self, sample_size):
    self.haps = {}                    # hap_seq : haplotype
    self.sample_size = sample_size
    self.freqs_initialized = False
    
  def __getitem__(self, key):
    return  self.haps[key]
    
  def __setitem__(self, key, value):
    self.haps[key] = value
    
  def __repr__(self):
    repr = ""
    for key in self.haps.keys():
      repr += '%s\t%s\n' % (key.__repr__(), self.haps[key].__repr__())
    
    totals = self.total_counts()
    repr += '\ntotal static count  %f   ' % totals[0]
    repr += 'total old count %f   ' % totals[1]
    repr += 'total new count %f   ' % totals[2] 
    return  repr
    
  def cmp_keys(self, left, right):
    left_name = self.haps[left].name
    right_name = self.haps[right].name
    result = 0
    if left_name < right_name:
      result = -1
    elif left_name > right_name:
      result = 1
      
    return  result

  def output(self):
    keys = self.haps.keys()
    keys.sort(self.cmp_keys)
    for key in keys:
      print '%s\t%f' % (self.haps[key].name,
                        self.haps[key].new_freq(self.sample_size)) 
    
  def has_key(self, key):
    return  self.haps.has_key(key) 
    
  # - Are new frequency estimates sufficiently close to old estimates
  #   for all haplotypes?
  #
  def converged(self, epsilon):
    conv = True
    for key in self.haps.keys():
      if not self.haps[key].converged(epsilon, self.sample_size):
        conv = False
        break
        
    return  conv
    
  def calc_freqs(self, phenos):
    self.reset()
    if self.freqs_initialized:           # Static counts have been set.
      for pheno_seq in phenos:
        if phenos[pheno_seq].ambiguous:
          for hap_set in phenos[pheno_seq].hap_sets:
          
            # - phenotype count * haplotype set weight
            #
            incr = phenos[pheno_seq].count * phenos[pheno_seq].hap_sets[hap_set]
            for hap in hap_set:
              self.haps[hap].iadd(incr)
          
    else:
      for pheno_seq in phenos:
        for hap_set in phenos[pheno_seq].hap_sets:
          incr = phenos[pheno_seq].count * phenos[pheno_seq].hap_sets[hap_set]
          for hap in hap_set:
            self.haps[hap].iadd(incr)
            if not phenos[pheno_seq].ambiguous:
              self.haps[hap].static_count = self.haps[hap].static_count + incr
          
      self.freqs_initialized = True
                
    
  # - Called at each iteration of EM algorithm.
  #
  def reset(self):
    for key in self.haps.keys():
      self.haps[key].reset()
      
  # - Called at start of EM algorithm.
  #
  def zero_freqs(self):
    for key in self.haps:
      self.haps[key].old_count = 0.0
      self.haps[key].new_count = 0.0
      
  # - Total haplotype frequencies.  Should always equal 1.0
  #
  def total(self):
    total = 0.0
    for key in self.haps.keys():
      total += self.haps[key].new_freq(self.sample_size)
      
    return  total
    
  # - Return total static, old and new counts.
  #
  def total_counts(self):
    total_static = 0.0
    total_old    = 0.0
    total_new    = 0.0
    for key in self.haps:
      total_static += self.haps[key].static_count
      total_old    += self.haps[key].old_count
      total_new    += self.haps[key].new_count
 
    return  total_static, total_old, total_new
    
    
    
if __name__ == '__main__':
  pass