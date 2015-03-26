#============================================================================
# File:      test_info.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 8/20/2                                             
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
# Copyright (c) 2002 R.C. Elston                                           
# All Rights Reserved                                                    
#============================================================================

import sagetest

class func_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def testparse(self):
    self.cmd = 'test_parse par > screen 2>&1'
    self.file_names = ['testparse.inf', 'testparse.out', 'screen']
    self.execute()

  def testfunc(self):
    self.cmd = 'test_func mypara myped1 mymld > out 2>&1'
    self.file_names = ['testfunc.inf', 'testfunc.out']
    self.execute()
    
  def testadj(self):
    self.cmd = '../pedinfo/pedinfo par ped > out 2>&1'
    self.file_names = ['pedinfo.inf']
    self.execute()
    
  def testbinary(self): 
    self.cmd = '../pedinfo/pedinfo par dat > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'screen']
    self.execute()
    
  def testtimer(self):
    self.cmd = 'testtimer.script > out 2>&1'
    self.file_names = ['timer.out', 'out']
    self.execute()

  def test_list(self):
    self.cmd = '../assoc/assoc par ped > out 2>&1 '
    self.file_names = ['assoc.inf', 'out']
    self.success_expected=0
    self.execute()
    
  def test_tai(self):
    self.cmd = '../freq/freq par ped > out 2>&1 '
    self.file_names=['freq.inf', 'out']
    self.execute()
    
  def test_poo(self): 
    self.cmd = '../pedinfo/pedinfo par dat > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'screen']
    self.execute()
    
    
    