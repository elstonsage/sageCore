#===================================================================
# File:      test_info.py
#
# Author:    Geoff Wedig
#
# History:   Updated from Old Version 2003-10-27
#
# Notes:     For use by test scripts (see src/c++/test_scripts)
#
# Copyright (c) 2003 R.C. Elston
# All Rights Reserved
#====================================================================

import sagetest
class parser_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.delta       = 0.001 # estimate difference between linux & windows is in .001.

  def test1(self):
    'sum, W3, means_test'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal params ped ibd >out 2>&1'
    self.delta       = 0.02
    self.execute()

  def test2(self):
    'W2'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal params ped ibd >out 2>&1'
    #self.delta       = 0.02
    self.execute()

  def test3(self):
    'default - cross product, mean_test'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal params ped ibd >out 2>&1'
    self.execute()

  def test4(self):
    'single vs. multiple, mean_test'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal params ped ibd >out 2>&1'
    self.delta       = 0.02
    self.execute()

  def test5(self):
    'sum, diff, W4, zero_marker, means_test'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'zero_test.treg', 'zero_test.treg_det', 'out']
    self.cmd         = 'sibpal g.par g.ped g.ibd >out 2>&1'
    self.execute()

  def test6(self):
    'DBH data - W3, means_test'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'abo.cor', 'abo.design', 'abo.qls', 'out']
    self.cmd         = 'sibpal params ped ibd >out 2>&1'
    self.execute()

  # ticket #700
  # Multiple IBD files.
  # The pair order is different in 3 ibd files.
  # 
  def test7(self):
    'Test for multiple IBD files'
    self.file_names  = ['sibpal.inf', 'noRVE.treg', 'noRVE.treg_det', 'out']
    self.cmd         = 'sibpal par ped ch16.ibd ch17.ibd ch18.ibd >out  2>&1'
    self.execute()

  def test8(self):
    'half sib test1 - W3, export_out, design & correlation matrix printout'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'traits.export', 'D3S1297.cor', 'D3S1297.design', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  def test9(self):
    'half sib test2 - W4, BLUP mean, export_out, data_dump, design & correlation matrix printout'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'traits.export', 'D3S1297.cor', 'D3S1297.design', 'D3S1297.qls', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  def test10(self):
    'mean_test - w1, export_out'
    self.file_names  = ['sibpal.inf', 'test.mean', 'test.mean_export', 'm1.mean', 'm2.mean', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  def test11(self):
    'mean_test, regression test with interaction'
    self.file_names  = ['sibpal.inf', 'means.out', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal param ped IBD1 >out 2>&1'
    self.execute()

  def test12(self):
    'pair-information file, data dump'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'D16S683.qls', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  def testX1(self):
    'x-linkage test'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'DXS6807.qls', 'DXS6807.design', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  def testX2(self):
    'x-linkage test with missing sex_info'
    self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'out']
    self.cmd         = 'sibpal par ped ibd >out 2>&1'
    self.execute()

  # Not included into daily build test due to time.
  # Testing is done whenever the code changes.
  #
  #def test_empp(self):
  #  'empirical p-value test'
  #  self.file_names  = ['sibpal.inf', 'traits.out', 'traits.det', 'out']
  #  self.cmd         = 'sibpal par1 ped1 ibd1 >out 2>&1'
  #  self.execute()
