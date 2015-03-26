#============================================================================
# File:      test.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 7/31/2                                             
#                                                                          
# Notes:     module for test scripts.
#                                                                          
# Copyright (c) 2002 R.C. Elston                                           
# All Rights Reserved                                                    
#============================================================================

import sys
import string

SUCCEEDED = 0
FAILED = 1
ABORTED = 2

TEST_CMD = 0
FILE_NAMES = 1
PURPOSE = 2
KEEP_FILES = 3

STRIP_LIST = ('S.A.G.E. v', 'COPYRIGHT', 'File generated on')

# - Information needed by a test script.
#
class info:
  def __init__(self, test_cmd, file_names, purpose):
    self.test_cmd = test_cmd
    self.file_names = file_names
    self.purpose = purpose

# - Display beginning of test message.
#
def begin(test_dir):
  print '\n'
  print '******************************************************************'
  print 'Running test, %s ...\n' % test_dir
  print 
  
# - Display message and exit.
#
def exit(exit_status, message, test_dir):
  print
  
  if exit_status == ABORTED:
    print '... test, %s, ABORTED.  ' % test_dir,
  elif exit_status == FAILED:
    print '... test, %s, FAILED.  ' % test_dir,
  elif exit_status == SUCCEEDED:
    print '... test, %s, SUCCEEDED.  ' % test_dir,
    
  print message
    
  print '******************************************************************'
  print '\n'

  sys.exit(exit_status)

# - Strip the indicated file of lines containing strings in strip_list.
#
def strip(file_name, strip_list):    
  try:
    file = open(file_name, 'r')
  except IOError:
    print 'Could not open file %s.' % file_name
    sys.exit(1)
    
  lines = file.readlines()
  file.close()
  
  for phrase in strip_list:
    lines = [line for line in lines if string.find(line, phrase) == -1]
 
  try:
    file = open(file_name + ".stripped", 'w+')
  except IOError:
    print 'Could not open file %s.' % file_name
    sys.exit(1)
    
  file.writelines(lines)
  file.close()

    
if __name__ == '__main__':
  strip(sys.argv[1], STRIP_LIST)      
    