#!/usr/local/bin/python

import sys
import os
import time
import stat
import string
import unittest
import common


seddir = os.path.abspath(sys.path[0])

class SAGE_Test (unittest.TestCase):
  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def  __init__(self, methodName="runTest", test_dir="",params=[], common_dir=''):
    unittest.TestCase.__init__(self, methodName)
    
    if test_dir == "":
      self.test_dir = methodName
    else:
      self.test_dir = test_dir

    if len(params):
      self.cmd        = params[0]
      self.file_names = params[1]

      self.doc = params[2]
    else:
      self.cmd        = ""
      self.file_names = []

      self.doc = unittest.TestCase.shortDescription(self)

    self.common_path   = common_dir
    self.relative_path = ""
    self.keep_files    = 0
    self.success_expected = 1

    self.delta   = 0.0001
    self.epsilon = 1.0e-10 

  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def shortDescription(self):
    doc = self.doc
    doc = doc and string.strip(string.split(doc, "\n")[0]) or None

    if doc != None:
      return self.test_dir + " : " + doc
    else:
      return None

  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def tearDown(self):
    os.chdir(self.execution_path)

  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def runTest(self):
    self.execute()
    
  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def execute(self):
    platform   = common.get_os()
    module     = os.path.split(os.path.abspath(""))[1]
    test       = self.test_dir
    start_time = time.strftime("%Y.%m.%d %H:%M:%S")

    return_value = self.perform_test();
    
    end_time = time.strftime("%Y.%m.%d %H:%M:%S")
    
    self.record_results(platform, module, test, start_time, end_time, return_value[0])
    
    self.failIf(return_value[0] != 0, return_value[1])

  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def perform_test(self):
    
    self.execution_path = os.path.abspath("")
    
    test_path = self.get_test_path()

    if (os.path.isdir(test_path) == 0):
      return (3, 'Could not find directory %s' % test_path)
    
    # Get the current directory and apply any relative pathing to it

    curr_path = os.path.abspath(self.relative_path)

    os.chdir(test_path)

    # - Execute test program.
    #


    test_return = os.system('%s/%s' % (curr_path, self.cmd))

    # - Test for success or failure appropriately
    #
    if (self.success_expected):
      if(test_return != 0):
        return (1, 'Test program execution failure.')
    else:
      if(test_return == 0):
        return (2, 'Test program execution non-failure.')
      
    # - Diff result outputs of test program against expected results.
    #
    redirect = '>'

    os.system('rm brokenfiles 2>/dev/null')

    for idx in range(len(self.file_names)):
      if idx != 0:
        redirect = '>>'

      file_name = self.file_names[idx]

      # Determine the expected file name - This is either the file name
      # with '.exp' appended or '.<OS>.exp' appended if there is an OS
      # specific version of the expected file

      OS = common.get_os()

      exp_file_name = '%s.%s.exp' % (file_name, OS)
      exp_file_path = './%s' % exp_file_name

      if os.access(exp_file_path, os.F_OK) == 0 :
        exp_file_name = '%s.exp' % file_name
        exp_file_path = './%s' % exp_file_name

      if(os.path.exists(exp_file_path) == 0):
        return (4, 'File %s not found.' % exp_file_name)

      # - Create a copy of the file as <file>.stripped that can be compared
      #   against the .exp file.  This is to preserve the existing file.
      #
      strip(file_name, STRIP_LIST)

      # - Make 'NANQ' and 'INF' uniform across platforms.
      #

      sed_run_cmd = 'sed -f %s/sed_script %s > %s/temp' % (seddir, file_name + ".stripped", seddir)
      os.system(sed_run_cmd)
      mv_cmd = 'mv %s/temp ./%s' % (seddir, file_name + ".stripped")
      os.system(mv_cmd)

      # do the diffing.

      command = "numdiff -d%s -e%s " % (self.delta, self.epsilon)

      ret_type = os.system('%s %s %s %s testdiffs' % \
            (command, exp_file_name, file_name + ".stripped", redirect))

      if ret_type == 0:
        # If there are no diffs, we can delete the .stripped file
        os.system('rm ' + file_name + '.stripped')

        if self.keep_files == 0:
          os.system('rm ' + file_name)
      else:
        os.system('echo ' + file_name + ' >>brokenfiles')


    results = os.stat('testdiffs')
    if(results[stat.ST_SIZE] != 0):
      return (5, 'Test Failure.  See %s/testdiffs for details.' % test_path)
    else:
      os.system('rm testdiffs')
      
    return (0, "Success")

  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def get_test_path(self):
    
    test_path = self.test_dir

    if self.common_path != '':
      test_path = '%s/%s' % (self.common_path, self.test_dir)

    return test_path
    
  ## ---------------------------------------------------------------------------
  ## Method
  ## ---------------------------------------------------------------------------
  def record_results(self, platform, module, test, start, end, result):
    try:
      timestamp = os.environ["SAGE_TIME_STAMP"]
      
      cmd = 'build_test_log db "' + timestamp + '" "' \
                                   + platform  + '" "' \
                                   + module    + '" "' \
                                   + test      + '" "' \
                                   + start     + '" "' \
                                   + end       + '"  ' \
                                   + str(result)
                                  

      os.system(cmd)
               
    except:
      no_op = "foo"
      
    try:
      statsfile = os.environ["SAGE_TEST_STATS"]
      
      cmd = 'build_test_log stats "' + statsfile + '" "' \
                                     + platform  + '" "' \
                                     + module    + '" "' \
                                     + test      + '" "' \
                                     + start     + '" "' \
                                     + end       + '"  ' \
                                     + str(result)
                                  
      os.system(cmd)
      
    except:
      no_op = "foo"
               
                                       

STRIP_LIST = ('S.A.G.E. v', 'COPYRIGHT', 'File generated on',
              'Remember you have agreed to add an appropriate statement (including the',
              'NIH grant number) under "acknowledgments" in any publication of results',
              'obtained by using this program. Suggested wording is:',
              '(Some of)The results of this paper were obtained by using the software',
              'package S.A.G.E., which was supported by a U.S. Public Health Service',
              'Resource Grant (RR03655) from the National Center for Research Resources.',
              '********************************************************************',
              '*****                                                          *****',
              '*****    Always check the INF file before viewing results.     *****'
             )


## ---------------------------------------------------------------------------
## Function: Strip the indicated file of lines containing strings in strip_list.
## ---------------------------------------------------------------------------
def strip(file_name, strip_list):    
  try:
    file = open(file_name, 'r')
  except IOError:
    print 'Could not open file %s for input.' % file_name
    sys.exit(1)
    
  lines = file.readlines()
  file.close()
  
  for phrase in strip_list:
    lines = [line for line in lines if string.find(line, phrase) == -1]
 
  try:
    file = open(file_name + ".stripped", 'w+')
  except IOError:
    print 'Could not open file %s for output.' % file_name
    sys.exit(1)
    
  file.writelines(lines)
  file.close()
