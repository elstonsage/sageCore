#============================================================================
# File:      locus.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   12/4/3  Created.
#            12/18/3 Modified to allow arbitrary number of alleles at a locus.  
#                                                                          
# Notes:     
#
#============================================================================

class genotype:
  def __init__(self, allele1, allele2):
    self.allele1 = allele1
    self.allele2 = allele2
    self.homozygous = (self.allele1 == self.allele2)
    self.heterozygous = not self.homozygous
    
  def __getitem__(self, index):
    allele = None
  
    if index == 0:
      allele = self.allele1
    elif index == 1:
      allele = self.allele2
      
    return  allele
    
  def __repr__(self):
    return  '[%d, %d]' % (self.allele1, self.allele2)

# - Maps allele names (strings) to indices and vice versa.
#
class locus:
  def __init__(self):
    self.names = {}
    self.indices = {}
    self.genos = []
    
  def __repr__(self):
    return  '%s\n%s\n%s' % (self.names.__repr__(), 
                            self.indices.__repr__(),
                            self.genos.__repr__()   )
    
  def add_allele(self, name):
    if not self.indices.has_key(name):
      new_index = len(self.indices)
      self.indices[name] = new_index
      self.names[new_index] = name
      
  def build_genotypes(self):
    count = self.a_count()
    for first in range(count):
      for second in range(first, count):
        self.genos.append(genotype(first, second))
      
  def a_name(self, index):
    return  self.names[index]
    
  def a_index(self, name):
    return  self.indices.get(name, None)
    
  def a_count(self):
    return  len(self.indices)
    
  def g_name(self, index):
    return  '%s/%s' % (self.a_name(self.genos[index][0]),
                       self.a_name(self.genos[index][1]) )
                       
  def g_alt_name(self, index):
    return  '%s/%s' % (self.a_name(self.genos[index][1]),
                       self.a_name(self.genos[index][0]) )  
                       
  def g_index(self, name):
    g_index = None
    for i in range(len(self.genos)):
      if self.g_name(i) == name or self .g_alt_name(i) == name:
        g_index = i
        break
        
    return  g_index
  
  def g_count(self):
    return  len(self.genos)
    

if __name__ == '__main__':
  l = locus()
  l.add_allele('a')
  l.add_allele('b')
  l.add_allele('c')
  l.add_allele('d')
  
  l.build_genotypes()
  for geno in l.genos:
    print '%d/%d' % (geno[0], geno[1])