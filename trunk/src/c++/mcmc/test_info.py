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

COMMON_DIRECTORY = 'tests'

TESTS = { 'test0'        : ('test_mcmc par ped mld gen > out',
                           ['out'],
                            'marker_likelihhods test'),
          'test1'        : ('test_mcmc par ped mld gen > out',
                           ['out'],
                            'marker_likelihhods test'),
        }


