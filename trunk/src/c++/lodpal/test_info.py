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

class lodpal_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.delta       = 0.01    # Sets delta for comparing numbers (as fraction of first number).
    self.epsilon     = 1e-5    # Sets the epsilon for comparing numbers close to zero as zero.

  def test1(self):
    'one-parameter model vs. two-parameter model'
    self.test_dir = 'test1'
    self.file_names  = ['one.out', 'one.lod',
                        'two.out', 'two.lod',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal par ped ibd >out 2>&1' 
    self.execute()
    
  def test2(self):
    'parent-of-origin test'
    self.test_dir = 'test2'
    self.file_names  = ['one.out', 'one.lod',
                        'two.out', 'two.lod',
                        'maternal_fixed.out', 'maternal_fixed.lod',
                        'paternal_fixed.out', 'paternal_fixed.lod',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal 11.par 11.ped 11mp.ibd.new 2>&1 >out'
    self.execute()

  def test3(self):
    'pair information file test'
    self.test_dir = 'test3'
    self.file_names  = ['alz21A.out', 'alz21C.out',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal 21.par 21.dat 21.ibd 2>&1 >out'
    self.execute()
    
  def testX2(self):
    'X-linkage with pair information file test'
    self.test_dir = 'testX2'
    self.file_names  = ['base.xln',  'base.lod',
                        'cov_g.xln', 'cov_g.lod',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal par ped ibd 2>&1 >out'
    self.execute()

  def testX3(self):
    'various X-linkage tests'
    self.test_dir = 'testX3'
    self.file_names  = ['base_equal.xln',  'base_equal.lod',
                        'base_not_equal.xln', 'base_not_equal.lod',
                        'ageonset_equal.xln', 'ageonset_equal.lod',
                        'ageonset_not_equal.xln', 'ageonset_not_equal.lod',
                        'l2_equal.xln', 'l2_equal.lod',
                        'l2_not_equal.xln', 'l2_not_equal.lod',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal par ped ibd 2>&1 >out'
    self.execute()

  def test(self):
    'DSP contrast option test'
    self.test_dir = 'test'
    self.file_names  = ['with_subset.out', 'with_subset.lod',
                        'without_subset.out',
                        'one_contrast.out', 'two_contrast.out',
                        'lodpal.inf', 'out'  ]
    self.cmd         = 'lodpal par ped ibd 2>&1 >out'
    self.execute()
