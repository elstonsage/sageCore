#!/usr/bin/env python

import string
import os
import sys
import UserDict
import math

class EnumeratedSet(UserDict.UserDict):

  def __init__(self):
    UserDict.UserDict.__init__(self)
    self._items = []
    self.count = 0

  def __setitem__(self, item, value):
    raise "Use insert, not __setitem__"

  def insert(self, item):
    if not self.data.has_key(item):
      self._items.append(item)
      self.data[item] = self.count
      self.count = self.count + 1

  def items():
    return self._items


def get_line(file):

  if not file:
    return None
  line = file.readline()

  while line and len(line) and line[0] == '#':
    line = file.readline()

  return line

def read_ibdfile(file_name, ibd_map, markers):

  ibd_file  = open(file_name, "r")

  file_markers = []

  line = string.strip(get_line(ibd_file))

  if line[:12] != "IBD File 1.9" and line[:12] != "IBD File 2.0" :
    print "Bad IBD file!  Header =",line
    sys.exit(1)

  line = string.strip(get_line(ibd_file))

  if line != "ANALYSIS":
    print "Bad IBD file!",line
    sys.exit(1)

  while 1:
    line = string.strip(get_line(ibd_file))
    if not len(line):
      break

  line = string.strip(get_line(ibd_file))

  if line != "MARKERS":
    print "Bad IBD file!",line
    sys.exit(1)

  while 1:
    line = string.strip(get_line(ibd_file))
    if not len(line):
      break
    file_markers.append(line)
    markers.insert(line)

  while 1:
    line = get_line(ibd_file)
    
    if not line or not len(line):
      break

    fields = string.split(line, ",")

    if len(fields) < 4:
      print "Bad field =", line
      continue

    ped = string.strip(fields[0])
    id1 = string.strip(fields[1])
    id2 = string.strip(fields[2])
    
    if id1 > id2:
      temp = id1
      id1  = id2
      id2  = temp

    map_id = "%s,%s,%s" % (ped, id1, id2)
    
    mapped_ibds = map(string.strip, string.split(string.strip(fields[3]), "  "))
    
    ibds = []
    for i in mapped_ibds :
      if i != "-----------------" :
        ibds = ibds + [float(i)]
      else :
        ibds = ibds + [-1]
        
    for i in range( min(len(ibds)/3, len(file_markers)) ):
      if not ibd_map.has_key(map_id):
        ibd_map[map_id] = {}

      ibd_map[map_id][ file_markers[i] ] = ( ibds[3*i], ibds[3*i+2] )


try:
  file1 = sys.argv[1]
  file2 = sys.argv[2]

except IndexError:
  print "Usage: [filenames]"
  sys.exit(1)

ibd_map1 = {}
ibd_map2 = {}
markers = EnumeratedSet()

read_ibdfile(file1, ibd_map1, markers)
read_ibdfile(file2, ibd_map2, markers)

se = 0.0
count = 0

for i in ibd_map1.keys():

  if not ibd_map2.has_key(i):
    continue

  for m in ibd_map1[i].items():

    mrk  = m[0]
    ibd1 = m[1]

    if not ibd_map2[i].has_key(mrk):
      continue

    ibd2 = ibd_map2[i][mrk]
    
    if ibd1[0] < 0 or ibd1[1] < 0 or ibd2[0] < 0 or ibd2[1] < 0 :
      continue

    er = (ibd1[0]-ibd2[0])**2 + (ibd1[1]-ibd2[1])**2 \
       + (ibd1[0]+ibd1[1]-ibd2[0]-ibd2[1])**2

    if er > .2:
      print "Error (%6.4f) %12s: %12s = (%6.4f,%6.4f), (%6.4f,%6.4f)" \
                              % (er, i,mrk,ibd1[0],ibd1[1],ibd2[0],ibd2[1])

    se = se + er

    count = count + 1

mse = se / count / 3

print mse
print count


