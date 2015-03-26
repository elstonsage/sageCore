#============================================================================
# File:      locus.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/15/4  Created.
#                    Added some functions when estimate1 was created.
#                                                                          
# Notes:     For pooled version for haplotype freq. estimating program.
#
#============================================================================

class locus:
  """ A genotype name for a 2 person pool might be 'a/a/b/c'.
      The corresponding sequence might be (2, 1, 1) indicating 2 of allele 0,
      1 of allele and 1 of allele 2 """

  def __init__(self, pool_size):
    self.a_names   = {}         # Allele names keyed by indices.
    self.a_indices = {}         # Allele indices keyed by names.
    self.g_names   = {}         # Genotype names keyed by sequence. 
    self.g_seqs    = {}         # Genotype sequences keyed by allele name.
    self.g_size = 2 * pool_size    # Number of alleles in a genotype.    

  def add_allele(self, name):
    if not name in self.a_indices:
      new_index = len(self.a_indices)
      self.a_indices[name] = new_index
      self.a_names[new_index] = name
      
  def add_geno(self, alleles):
    seq = self.g_seq(alleles)
    if not seq in self.g_names:
      name = self.g_create_name(seq)
      self.g_names[seq] = name
      self.g_seqs[name] = seq
      
  def a_name(self, index):
    return  self.a_names[index]
    
  def a_index(self, name):
    return  self.a_indices[name]
    
  def a_count(self):
    return  len(self.a_names)
    
  def g_name(self, seq):
    return  self.g_names[seq]
    
  def g_seq(self, alleles):
    """ Count alleles by building a temporary dictionary whose keys are 
        allele names and whose values are lists whose lengths equal the 
        number of copies of the allele. """
    if alleles[0] == '?':
      return  None
        
    a_counts = {}
    for a in alleles:
      a_counts.setdefault(a, []).append(None)
      
    seq = []
    for index in range(self.a_count()):
      seq.append(len(a_counts.get(self.a_name(index), [])))

    return  tuple(seq)
    
  # - Return a list of allele indices.
  #
  def alpha_to_indices(self, alleles):
    index_list = []
    for allele in alleles:
      index_list.append(self.a_index(allele))
      
    return  tuple(index_list)          
    
  def g_create_name(self, g_seq):
    l_name = []
    for index in range(len(g_seq)):
      for a in range(g_seq[index]):
        l_name.append(self.a_name(index)) 
        
    l_name.sort()
    str_name = ""
    for a in l_name:
      str_name += a + '/'
      
    return  str_name[:-1]
    
  def g_alt_create_name(self, allele_seq):
    name = ''
    for a_index in allele_seq:
      name += '%s/' % self.a_name(a_index)
      
    return  name[:-1]
    

if __name__ == '__main__':
  loc = locus(2)
  for a in ('d', 'b', 'c', 'a'):
    loc.add_allele(a)
    
  loc.add_geno(['a', 'a', 'a', 'b'])
  loc.add_geno(['b', 'b', 'b', 'a'])
  
  for name in loc.g_seqs:
    print name
    
  print loc.alpha_to_indices(['a', 'b'])
