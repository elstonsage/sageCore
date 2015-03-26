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

class fcor_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  def test1(self):
    'subtypes'
    self.test_dir = 'test1'
    self.file_names  = ['fcor.out', 'fcor.inf', 'out']
    self.cmd         = 'fcor data.par data.ped 2>&1 >out'
    self.execute()                
    
  def test2(self):
    'sub & main types, pair files'
    self.test_dir = 'test2'
    self.file_names  = ['fcor.out',
                        'fcor.pair', 'fcor.inf', 'out' ]
    self.cmd         = 'fcor data2.par data2.ped 2>&1 >out'
    self.execute()

  def test3(self):
    'homogeneity test, variance-covariance, alt file'
    self.test_dir = 'test3'
    self.file_names  = ['fcor.cov',   'fcor.out',
                        'fcor.inf', 'out' ]
    self.cmd         = 'fcor data3.par data3.ped 2>&1 >out'
    self.execute()                        
    
  def test4(self):
    'simulated data'
    self.test_dir = 'test4'
    self.file_names  = ['TEST1.cov', 'TEST1.out', 'TEST1.det',
                        'TEST1.alt', 'TEST1.pair',
                        'fcor.inf', 'out' ]
    self.cmd         = 'fcor simul.par simul.dat 2>&1 >out'
    self.execute()                    
    
  def test7(self):
    'with locus file'
    self.test_dir = 'test7'
    self.file_names  = ['fcor.out', 'fcor.inf',
                        'genome.inf', 'out' ]
    self.cmd         = 'fcor params ped loc 2>&1 > out'
    self.delta       = 0.02
    self.execute()

