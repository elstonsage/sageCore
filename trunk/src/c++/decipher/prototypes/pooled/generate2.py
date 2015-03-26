#============================================================================
# File:      generate2.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   11/18/4.  Created.
#                                                                          
# Notes:     generate haplotype combinations consistent with phenotypes
#            for pool data.  Haplotypes are represented by a tuple of
#            positive integers.  phenotypes are represented by a sequence of
#            sequences each inner sequence represents one locus.  It's ele-
#            ments are allele counts. 
#
#============================================================================

from copy import copy, deepcopy

def generate(pheno_seq, hap_seq_sets):
  build(init_comb(pheno_seq[0]), pheno_seq, hap_seq_sets)

def build(comb, pheno_seq, hap_seq_sets):
  if(len(comb[-1]) == len(pheno_seq)):
    hap_seq_sets.append(tuplize(comb))
    
  else:
    start = find_start(comb)
    locus = find_current_locus(comb, start)
    uc = uncovered_count(comb, len(pheno_seq))
    aa = available_alleles(comb, pheno_seq, start, locus)
  
    allele_groups = []
    generate_sets(aa, [], uc, allele_groups)
    
    for group in allele_groups:
      new_comb = deepcopy(comb)
      for a in range(len(group)):
        new_comb[start + a].append(group[a])
        
      build(new_comb, pheno_seq, hap_seq_sets)
      
# - Return list whose elements are tuples.
#
def tuplize(comb):
  new_comb = []
  for elem in comb:
    new_comb.append(tuple(elem))
    
  return  new_comb
  
  
def init_comb(pheno):
  comb = []
  for l in range(len(pheno)):
    for rep in range(pheno[l]):
      comb.append([l])
      
  return  comb
  
#   0011 yeilds 2,   0001 yeilds 1,   etc.
#                    11
#  
# if combination is complete, None is returned.
#
def uncovered_count(comb, locus_count):
  if len(comb[-1]) == locus_count:
    return  None
  else:
    start = find_start(comb)
    return  comb[start:].count(comb[start])
  
  
# - Return index of first sublist that is shorter than the first one.
#
def find_start(comb):
  start = 0
  length = len(comb[0])
  for h in range(len(comb)):
    if len(comb[h]) < length:
      start = h
      break
      
  return  start
  
# - Return index of next locus if start is 0 or current locus if start is
#   not 0.
#
def find_current_locus(comb, start):
    return  len(comb[start])
  

# - Enumerate alleles at the given locus.
#
def locus_alleles(pheno):
  alleles = []
  for l in range(len(pheno)):
    for rep in range(pheno[l]):
      alleles.append(l)
      
  return  alleles
  
# - Return 'unused' alleles in the current locus.
#
def available_alleles(comb, pheno_seq, start, current_locus):
  alleles = locus_alleles(pheno_seq[current_locus])
  for hap in comb[:start]:
    alleles.remove(hap[-1])
    
  return  alleles
  
# - find all unique, unordered groups of size, set_size, starting with
#   root, utilizing elements in pool without replacement.
#
def generate_sets(pool, root, set_size, results):
  if len(root) < set_size:
    for elem in pool:
      new_pool = copy(pool)
      new_pool.remove(elem)
      new_root = copy(root)
      new_root.append(elem)
      generate_sets(new_pool, new_root, set_size, results)
  else:
    root.sort()
    if not root in results:
      results.append(root)
 

    
if __name__ == '__main__':
  pheno_seq = [(2, 2), (2, 2), (1, 3)]
  hap_seq_sets = []
  
  generate(pheno_seq, hap_seq_sets)
  
  for set in hap_seq_sets:
    print set
  
