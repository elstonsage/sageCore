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

class LSF_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test(self):
    'LSF Testing'
    self.cmd = 'LSFtest >LSFtest.out'
    self.file_names = ['LSFtest.out' ]
    self.execute()

