#============================================================================
# File:      test_info.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 7/31/2                                             
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


class lvec_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_ivg(self):
    'Tests of Meiosis maps, inheritance vector generation and nodes'
    self.cmd = 'test_ivg params ped loc 2>&1 >ivg.out'
    self.file_names = [ 'ivg.out' ]
    self.execute()

  def test_nodes(self):
    'Tests of the nodes'
    self.cmd = 'test_nodes >test_nodes.out'
    self.file_names = [ 'test_nodes.out' ]
    self.execute()

  def test_fbc(self):
    'Tests of the fbc'
    self.cmd = 'test_ivg params ped loc 2>&1 >out'
    self.file_names = [ 'out' ]
    self.execute()
