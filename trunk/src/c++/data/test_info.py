#============================================================================
# File:      test_info.py
#
# Author:    Geoff Wedig
#
# History:   Updated from Old Version 2003-10-27
#
# Notes:     For use by test scripts (see src/c++/test_scripts).
#
# Copyright (c) 2003 R.C. Elston
# All Rights Reserved
#============================================================================

import sagetest

class data_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  def test_ruleset(self):
    """
    Purpose:  test ArgumentRuleset and attendant classes.
    
    Basis:    hand inspection. 
    """
    
    self.cmd = 'test_cmdline > screen 2>&1'
    self.file_names = ['screen']
    self.execute()
    
