#============================================================================
# File:      tree.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/21/4.  Created.
#            1/22/4.  Set aside as ill conceived.
#                                                                          
# Notes:     Create and traverse a tree for enumerating all of possible
#            haplotype combinations given pool haplotype count and count of 
#            possible haplotypes. 
#
#============================================================================

from copy import copy

g_size = 4        # Number of haplotypes in a pool.
h_count = 4       # Number of possible haplotypes.


class node:

  def __init__(self, h_index, count, path_length):
    self.h_index = h_index          # Haplotype index.  Tree depth.
    self.count = count              # Count of haplotype represented by h_index.
    self.path_length = path_length  # Number of haplotypes represented by each path.
    print self
    
    if self.h_index < h_count - 1:
      self.left_child = node(self.h_index + 1, 0, self.path_length)
    else:
      self.left_child = None
      
    if path_length < g_size:
      self.sibling = node(self.h_index, self.count + 1, self.path_length + 1)
    else:
      self.sibling = None
      
  def __repr__(self):
    return  '%d/%d/%d' % (self.h_index, self.count, self.path_length)
      

def hap_sets(node, sets, set):
  """ Traverse tree from node to leaves building a list of lists where
      each index of inner list represent a haplotype and values represent a 
      count of that haplotype """
      
  if node.left_child:
    if node.sibling:
      hap_sets(node.left_child, sets, set)
    else:
      set.append(node.count)
      hap_sets(node.left_child, sets, copy(set))
      
  if node.sibling:
    hap_sets(node.sibling, sets, set)
  else:
    set.append(node.count)
    sets.append(set)
    
    
    
if __name__ == '__main__':
  root = node(0, 0, 0)

  sets = []
  hap_sets(root, sets, [])
  
#  for set in sets:
#    print set  