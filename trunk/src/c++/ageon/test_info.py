## -----------------------------------------------------------------------------
## Begin File: test_info.py
## -----------------------------------------------------------------------------
## -------------------------------------------------------------------------- ##
##                                   AGEON                                    ##
##                                Test Script                                 ##
## -------------------------------------------------------------------------- ##
# Author:    Kevin Cartier
#
# History:   Last revised on 23 May 2007
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
class ageon_tests (sagetest.SAGE_Test):

   def setUp(self):
     self.common_path = 'tests'


   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test1(self):
      'standard test'
      self.cmd        = 'ageon -p ageon.par -d ageon.ped >out 2>&1'
      self.file_names = ['1.det', '1.sum', '1.ped', '1.par', 'out']
      self.epsilon    = .0001
      self.delta      = 0.01
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test2(self):
      'multiple analyses test'
      self.cmd        = 'ageon par ped >out 2>&1'
      self.file_names = ['ageon_analysis1.det', 'ageon_analysis1.sum',
                         'ageon_analysis2.det', 'ageon_analysis2.sum',
                         'ageon_analysis3.det', 'ageon_analysis3.sum',
                         'ageon_analysis4.det', 'ageon_analysis4.sum',
                         'ageon.inf', 'out']
      self.epsilon    = .0001
      self.delta      = 0.01
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test_audrey(self):
      'one-column vs. two-column test with dbh data'
      self.cmd        = 'ageon par ped >out 2>&1'
      self.file_names = ['ageon_analysis1.det', 'ageon_analysis1.sum',
                         'ageon_analysis2.det', 'ageon_analysis2.sum',
                         'ageon.inf', 'out']
      self.epsilon    = .0001
      self.delta      = 0.01
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test_extra_output(self):
      self.cmd = 'ageon -p ageon.par -d ageon.ped >/dev/null 2>&1; ../../../app/test_app new.par ageon.ped >/dev/null 2>&1'
      self.file_names = ['ageon.inf']
      self.execute()

   ## --------------------------------------------------------------------------
   ## Method
   ## --------------------------------------------------------------------------
   def test_pool(self):
      'pooling test'
      self.cmd        = 'ageon -p ageon.par -d ageon.ped >out 2>&1'
      self.file_names = ['ageon_analysis1.det', 'ageon_analysis1.sum',
                         'ageon_analysis1.ped', 'ageon_analysis1.par',
                         'ageon.inf', 'out']
      self.epsilon    = .0001
      self.delta      = 0.01
      self.execute()


## -----------------------------------------------------------------------------
## End File: test_info.py
## -----------------------------------------------------------------------------
