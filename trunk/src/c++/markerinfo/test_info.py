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

class markerinfo_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['markerinfo.out', 'out',
                        'markerinfo_clean.ped', 'markerinfo_clean.par',
                        'markerinfo.inf', 'genome.inf'  ]
    self.cmd         = 'markerinfo -p param -d ped -l loc >out 2>&1'

  def test1(self):
    'misc tests'
    self.test_dir = 'test1'
    self.execute()
    
  def test2(self):
    'misc tests'
    self.test_dir = 'test2'
    self.execute()

  def test3(self):
    'misc tests'
    self.test_dir = 'test3'
    self.execute()
    
  def test4(self):
    'misc tests'
    self.test_dir = 'test4'
    self.file_names  = ['markerinfo.out', 'out',
                        'markerinfo.inf', 'genome.inf'  ]
    self.execute()

  #def test5(self):
  #  'misc tests'
  #  self.test_dir = 'test5'
  #  self.cmd         = 'markerinfo -p par -d ped -l loc >out 2>&1'
  #  self.execute()

  def testX(self):
    'X chromosome test'
    self.test_dir = 'testX'
    self.execute()

  def testY(self):
    'Y chromosome test'
    self.test_dir = 'testY'
    self.execute()

  def test_sibs(self):
    'sib pairs test'
    self.test_dir = 'test_sibs'
    self.file_names  = ['markerinfo.out', 'out',
                        'markerinfo.inf', 'genome.inf'  ]
    self.execute()

  def test_noncodom(self):
    'non-codominant marker testing'
    self.test_dir = 'test_noncodom'
    self.file_names  = ['markerinfo.out', 'out',
                        'markerinfo.inf', 'genome.inf'  ]
    self.execute()
