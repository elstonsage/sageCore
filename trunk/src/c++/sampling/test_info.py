#============================================================================
#
#	THIS COPY SPECIALIZED FOR ASSOC by Stephen Gross 10 Apr 2003
#
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

COMMON_DIRECTORY = ''

import sagetest

class sampling_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_sampling(self):
    'Basic sampling test'
    self.cmd = 'test_sampling par.txt dbhcomt.dat > out'
    self.file_names = ['out']
    self.execute()

