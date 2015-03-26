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

class maxfunapi_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test2(self):
    'Test using MAXFUN::Function'
    self.cmd = "maxtest2 > out"
    self.file_names = ['out']
    self.epsilon = 0.0001
    self.delta = 0.01
    self.execute()
