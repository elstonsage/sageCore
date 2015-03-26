#============================================================================
# File:      combinations.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/20/4.  Created.
#                                                                          
# Notes:     Generate all combinations by choosing one item from each set in
#            the sequence. 
#
#============================================================================

def combinations(set_seq):
  nested_combs = reduce(lambda s1, s2: [(x,y) for x in s1 for y in s2], set_seq)
  combs = []
  for item in nested_combs:
    combs.append(tuple(flatten(item)))
    
  return  combs
  
def flatten(nested, result = None):
  """ See 'Python Cookbook', 2002.  P 23. """
  
  if result is None:
    result = []
    
  for item in nested:
    if type(item) is int:
      result.append(item)
    else:
      flatten(item, result)
      
  return  result
  

if __name__ == '__main__':
  set_seq = ((1, 2), (4, 5, 6), (7, 8))
  
  for item in combinations(set_seq):
    print item