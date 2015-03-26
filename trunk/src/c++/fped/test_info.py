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

class rped_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_filter(self):
    'Multipedigree Filter Tests'
    self.cmd = "filter_test filter.dat '(2(A4,1X),T21,A1,T10,2(1X,A4),T24,A1)' > out"
    self.file_names = ['out']
    self.execute()
