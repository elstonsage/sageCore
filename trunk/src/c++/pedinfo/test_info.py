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

class pedinfo_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test1(self):
    self.cmd = 'pedinfo dbh.par dbh.dat > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'pedinfo.out', 'screen']
    self.execute()

  def test2(self):
    self.cmd = 'pedinfo mypara myped1 > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'analysis1.out', 'analysis2.out', 'screen']
    self.execute()
    
  def test_loops(self):
    'non-marriage loop counting'
    self.cmd = 'pedinfo par ped > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'pedinfo.out', 'screen']
    self.execute()
    
  def test_cons_loops(self): 
    'find consanquineous mating pairs'
    self.cmd = 'pedinfo par ped > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'pedinfo.out', 'screen']
    self.execute()
    
  def test3(self):
    'added to test new "base_trait" information'
    self.cmd = 'pedinfo par ped > screen 2>&1'
    self.file_names = ['pedinfo.inf', 'a1.out', 'a2.out', 'a3.out', 'screen']
    self.execute()
