#============================================================================
# File:      count.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/20/4.  Created.
#                                                                          
# Notes:     Generate successive numbers using an arbitrary base. 
#
#============================================================================


def count(places, base):
  """ Note: least significant digit is the leftmost. """
  
  initial_call = True
  number = [0] * places
  
  if initial_call:
    initial_call = False
    yield  tuple(number)
    
  while number != [base - 1] * places:
    incr(number, base - 1)
    yield  tuple(number)
    
    
def incr(number, max_digit):
  for place in range(len(number)):
    while number[place] != max_digit:
      number[place] += 1
      return
    for lesser_place in range(place + 1):
      number[lesser_place] = 0

    
if __name__ == '__main__':
  for n in count(8, 2):
    print n