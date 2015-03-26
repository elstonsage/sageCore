## -----------------------------------------------------------------------------
## Begin File: test_info.py
## -----------------------------------------------------------------------------
## -------------------------------------------------------------------------- ##
##                                   TDTEX                                    ##
##                                Test Script                                 ##
## -------------------------------------------------------------------------- ##
# Author:    Kevin Cartier
#
# History:   Last revised on 22 May 2007
#
# Notes:     For use by test scripts (see src/c++/test_scripts).
#
#            While in parent of the test directory, use runall to run all tests.
#
#            While in parent of the test directory, use runtest <test directory>
#            to run a specific test.
#
#            While in test directory, use newexps to make the current
#            test output files the expected files.
#
# Copyright (c) 2007 R.C. Elston
# All Rights Reserved
## -------------------------------------------------------------------------- ##


## ************************************************************************** ##
## **************              Imported Modules                ************** ##
## ************************************************************************** ##
import sagetest


## ************************************************************************** ##
## **************              Class Definition                ************** ##
## ************************************************************************** ##
class tdtex_tests (sagetest.SAGE_Test) :

   def setUp(self):
      self.common_path = 'tests'

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test1(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test1'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test2(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test2'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test3(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test3'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test4(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test4'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test5(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test5'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test6(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test6'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test9(self):
      self.cmd = 'tdtex -p params -d ped  >display 2>&1'
      self.test_dir = 'test9'
      self.file_names = ['tdtex.inf', 'tdtex.out', 'display']
      self.execute()


## -----------------------------------------------------------------------------
## End File: test_info.py
## -----------------------------------------------------------------------------
