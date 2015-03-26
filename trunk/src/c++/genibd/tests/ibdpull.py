#!/usr/bin/env python

import string
import os
import sys
import UserDict
import math

def get_line(file):

  if not file:
    return None
  line = file.readline()

  while line and len(line) and line[0] == '#':
    line = file.readline()

  return line

def read_ibdfile(file_name, ibd_map):

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

    map_id = "%s\t%s\t%s" % (ped, id1, id2)
    
    mapped_ibds = map(string.strip, string.split(string.strip(fields[3]), "  "))
    
    ibds = []
    for i in mapped_ibds :
      if i != "-----------------" :
        ibds = ibds + [float(i)]
      else :
        ibds = ibds + [-1]
        
    for i in range( min(len(ibds)/3, len(file_markers)) ):
      if not ibd_map.has_key(file_markers[i]):
        ibd_map[file_markers[i]] = {}

      ibd_map[ file_markers[i] ][map_id] = ( ibds[3*i], ibds[3*i+2] )


try:
  file1 = sys.argv[1]
  marker = sys.argv[2]

except IndexError:
  print "Usage: [filenames]"
  sys.exit(1)

ibd_map1 = {}

read_ibdfile(file1, ibd_map1)

for i in ibd_map1[marker]:
  val = "%s\t%f\t%f" % (i, ibd_map1[marker][i][0], ibd_map1[marker][i][1]) 
  print val
  

