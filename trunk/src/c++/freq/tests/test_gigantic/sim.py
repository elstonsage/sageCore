#!/usr/local/bin/python

import random

print "id",

for i in range(0, 1000):
  print "marker" + str(i),

for i in range(0, 1000):
  print str(i),
  
  

  for j in range(0, 1000):
  
    marker = ""
    
    if(random.randint(0, 1) == 0):
      marker += "A"
    else:
      marker += "B"
      
    marker += "/"
    
    if(random.randint(0, 1) == 0):
      marker += "A"
    else:
      marker += "B"
    
    print marker,
    
  print
  