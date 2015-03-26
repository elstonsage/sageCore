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

# - test5 intentionally omitted.  Must use local runtest5 script for this test.
# - test10 also omitted

class ped_calc_tests (sagetest.SAGE_Test) :
  def setUp(self):
    self.common_path = 'tests'
    self.file_names = ['out']

  def test_bin_pen(self):
    'Tests BinaryPenetranceCalculator'
    self.cmd = 'test_bin_pen_calc 2>&1 >out'
    self.execute()

  def test_fam_resid_adj(self):
    'Tests Family Residual Adjustment classes'
    self.cmd = 'test_fam_resid_adj 2>&1 >out'
    self.execute()

