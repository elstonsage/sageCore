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

class functors_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test1(self):
    'Basic test of functors'
    self.cmd = "test_functor > out"
    self.file_names = ['out']
    self.execute()
