#============================================================================
# File:      equation8.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 11/21/02
#                                                                          
# Notes:     Evaluates equation 8 of Cleves/Elston. 1997.
#
# Copyright (c) 2002 R.C. Elston                                           
# All Rights Reserved
#============================================================================

from math import log10

# - f  -> female recombination fraction.
#   m  -> male recombination fraction.
#   2n -> number of families.
#   a  -> number of female informative families of type I.
#   c  -> number of male informative families of type I.
#
def eq8(f, m, n, a, c):
  term1 = ((1.0 - f) ** 2.0 + f ** 2.0) ** a
  term2 = (2.0 * f * (1.0 - f)) ** (n - a)
  term3 = ((1.0 - m) ** 2.0 + m ** 2.0) ** c
  term4 = (2.0 * m * (1.0 - m)) ** (n - c)
  
  return term1 * term2 * term3 * term4
  
  
def lod(f, m, n, a, c):
  return  log10(eq8(f, m, n, a, c) / eq8(.5, .5, n, a, c))
  
def class1_like(theta):
  return  (1 - theta) * (1 - theta) + theta * theta
  
def class2_like(theta):
  return  2 * theta * (1 -theta)

def class1_lod(theta):
  return  log10(class1_like(theta) / class1_like(.5))
  
def class2_lod(theta):
  return  log10(class2_like(theta) / class2_like(.5))  
  