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

class mfsubmodels_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_type_specific_submodel(self):
    'TypeSpecificSubmodel test'
    self.cmd = 'test_type_specific_submodel par ped > out'
    self.file_names = ['out']
    self.epsilon=1e-5
    self.execute()

