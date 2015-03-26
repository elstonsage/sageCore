#!/usr/local/bin/python
#============================================================================
# File:      estimate
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/5/4.  Created.
#                                                                          
# Notes:     Parses data for haplotype estimating program.
#
#============================================================================

def parse(data):
  records = []
  lines = data.readlines()
  for line in lines[1:]:
    fields = line.split()
    records.append([ fields[0] ])
    for index in range(1, len(fields)):
      records[-1].append(fields[index].split('/'))
      
  return  records
    
    
if __name__ == '__main__':
  data = open('../tests/test1/sage.dat', 'r')
  records = parse(data)
  for record in records:
    print record