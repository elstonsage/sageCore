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

class reltest_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  def test_aneurysm(self):
    'misc tests'
    self.test_dir = 'test_aneurysm'
    self.file_names  = ['reltest.sum', 'out',
                        'reltest.inf', 'genome.inf'  ]
    self.cmd         = 'reltest params ped markers genome 2>&1 >out'
    self.execute()
    
  def test_audrey(self):
    'misc tests'
    self.test_dir = 'test_audrey'
    self.file_names  = ['reltest.sum', 'out',
                        'reltest.inf', 'genome.inf'  ]
    self.cmd         = 'reltest params ped locus genome 2>&1 >out'
    self.execute()

  def test_german(self):
    'misc tests'
    self.test_dir = 'test_german'
    self.file_names  = ['test.sum', 'out',
                        'reltest.inf', 'genome.inf'  ]
    self.cmd         = 'reltest params ped loc genome 2>&1 >out'
    self.delta       = 0.01
    self.execute()
    
  def test_lupus(self):
    'misc tests'
    self.test_dir = 'test_lupus'
    self.file_names  = ['reltest.sum', 'out',
                        'reltest.inf', 'genome.inf'  ]
    self.cmd         = 'reltest params ped markers genome 2>&1 >out'
    self.execute()
    
  #def test_snp(self):
  #  'test for snp data with a few pair'
  #  self.test_dir = 'test_snp'
  #  self.file_names  = ['reltest.sum', 'out',
  #                      'reltest.inf', 'genome.inf'  ]
  #  self.cmd         = 'reltest par chr01.ped loc gen 2>&1 >out'
  #  self.execute()
