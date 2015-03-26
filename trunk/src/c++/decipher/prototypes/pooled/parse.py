#============================================================================
# File:      parse.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   1/30/4.  Created.
#                                                                          
# Notes:     Parses POOLED data for haplotype estimating program.
#
#============================================================================

def parse(data):
  records = []
  lines = data.readlines()
  for line in lines[1:]:                  # Ignore header.
    fields = line.split()
    records.append([fields[0]])           # Individual id.
    for index in range(1, len(fields)):
      records[-1].append(fields[index].split('/'))
      
  return  records
    
    
if __name__ == '__main__':
  data = open('../tests/test5/sage.dat', 'r')
  records = parse(data)
  for record in records:
    print record