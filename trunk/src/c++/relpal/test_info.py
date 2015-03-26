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

class relpal_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.delta       = 0.01    # Sets delta for comparing numbers (as fraction of first number).
    self.epsilon     = 1e-5    # Sets the epsilon for comparing numbers close to zero as zero.

  def test1(self):
    'no genetic effect'
    self.test_dir = 'test1'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.export',
                        'relpal.inf', 'out'  ]
    self.cmd         = 'relpal par model1rep1.dat ibd11 >out 2>&1' 
    self.execute()
    
  def test2(self):
    'with genetic effect'
    self.test_dir = 'test2'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.export',
                        'relpal.inf', 'out'  ]
    self.cmd         = 'relpal par model2rep1.dat ibd21 >out 2>&1' 
    self.execute()

  def test3(self):
    'different pedigree structure - too little ped count'
    self.test_dir = 'test3'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.inf',
                        '1_loc2.dat', '1_loc2.debug', 'out'  ]
    self.cmd         = 'relpal par ped ibd_5types >out 2>&1' 
    self.execute()

  def test4(self):
    'multiple locations'
    self.test_dir = 'test4'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.inf',
                        'out'  ]
    self.cmd         = 'relpal par ped ibd >out 2>&1' 
    self.execute()

  def test5(self):
    'ibd_state file use, residual file'
    self.test_dir = 'test5'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.inf',
                        'relpal.resid', 'out'  ]
    self.cmd         = 'relpal par ped ibd >out 2>&1' 
    self.execute()

  def test_binary(self):
    'binary traits, first_level test'
    self.test_dir = 'test_binary'
    self.file_names  = ['relpal.out', 'relpal.det', 'relpal.export',
                        'relpal.inf', 'out'  ]
    self.cmd         = 'relpal par ped >out 2>&1'
    self.execute()
    
