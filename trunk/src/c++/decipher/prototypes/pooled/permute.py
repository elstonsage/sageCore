#============================================================================
# File:      permute.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   10/21/4.  Created.
#                                                                          
# Notes:     Permute list of elements which may or may not be unique. 
#
#============================================================================

from copy import copy

def permute(elements, start, results):
  """
  Exiter algorithm from www.bearcave.com/random_hacks/permute.html
  """
  
  if start == len(elements) - 1:
    if elements not in results:
      results.append(copy(elements))
  else:
    for i in range(start, len(elements)):
      elements[i], elements[start] = elements[start], elements[i]
      permute(elements, start + 1, results)
      elements[i], elements[start] = elements[start], elements[i]
  
    
if __name__ == '__main__':
  results = []
  
  permute([0,1,0,1], 0, results)
  
  for permutation in results:
    print permutation
