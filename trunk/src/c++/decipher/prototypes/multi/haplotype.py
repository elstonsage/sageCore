#============================================================================
# File:      haplotype.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   12/4/3  Created.
#            12/18/3 Modified to allow arbitrary number of allele at a locus.
#                                                                          
# Notes:     
#
#============================================================================

from locus import *

# - Stores haplotype definitions in allele sequence and string formats
#   as well as haplotype counts.  Constructor takes a sequence of 
#   indices indicating which allele is present at each locus.
#
class haplotype:

  # - Generate string representation of a haplotype.
  #
  def generate_name(self, hap_seq, loci):
    str = ""
    for l in range(len(hap_seq)):
      str += loci[l].a_name(hap_seq[l])
      
    return  str    

  def __init__(self, hap_seq, loci):
    self.name = self.generate_name(hap_seq, loci)
    self.old_count = 0.0
    self.new_count = 0.0
    self.static_count = 0.0   # - Portion of count not changing during 
                              #   algorithm iteration.
    
  def __repr__(self):
    return  '%s\tstatic count   %f\told count   %f\tnew count   %f' % (self.name, self.static_count, 
                                                                       self.old_count, self.new_count)
                                
  def iadd(self, value):
    self.new_count += value
    
  def reset(self):
    self.old_count = self.new_count
    self.new_count = self.static_count
    
    
  # - Are new counts sufficiently close to old counts?
  #
  def converged(self, epsilon, sample_size):
    return  abs(self.new_freq(sample_size) - self.old_freq(sample_size)) < epsilon
    
  def new_freq(self, sample_size):    
    return  self.new_count / sample_size / 2
    
  def old_freq(self, sample_size):    
    return  self.old_count / sample_size / 2

 
class haplotype_set:
  def __init__(self, sample_size):
    self.haps = {}
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
    if self.freqs_initialized:
      for p_key in phenos.keys():
        if phenos[p_key].ambiguous:
          for hp_key in phenos[p_key].hap_pairs:
            incr = phenos[p_key].count * phenos[p_key].hap_pairs[hp_key]
            self.haps[hp_key[0]].iadd(incr)
            self.haps[hp_key[1]].iadd(incr)
          
    else:
      for p_key in phenos.keys():
        for hp_key in phenos[p_key].hap_pairs:
          incr = phenos[p_key].count * phenos[p_key].hap_pairs[hp_key]
          self.haps[hp_key[0]].iadd(incr)
          self.haps[hp_key[1]].iadd(incr)
          if not phenos[p_key].ambiguous:
            self.haps[hp_key[0]].static_count = self.haps[hp_key[0]].static_count + incr
            self.haps[hp_key[1]].static_count = self.haps[hp_key[1]].static_count + incr
          
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
    
    
# - Convert from an integer representation of a haplotype phases
#   to a sequence representation.
#
def int_to_seq(i, length):
  seq = []
  seq.insert(0, i % 2)
  i /= 2
  while i != 0:
    seq.insert(0, i % 2)
    i /= 2
    
  # - Pad w. 0's if necessary.
  #
  while len(seq) < length:
    seq.insert(0, 0)
    
  return  seq    
        
    
if __name__ == '__main__':
  print int_to_seq(15, 4)    