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

class pairs_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  def test1(self):
    'relmatrix test1'
    self.test_dir = 'test1'
    self.file_names  = ['test.out', 'test.inf', 'out'  ]
    self.cmd         = 'testrelmatrix data3.par data3.ped >out 2>&1' 
    self.execute()
    
  def test2(self):
    'relmatrix test2'
    self.test_dir = 'test2'
    self.file_names  = ['test.out', 'test.inf', 'out'  ]
    self.cmd         = 'testrelmatrix data2.par data2.ped >out 2>&1' 
    self.execute()

  def test3(self):
    'relmatrix test with unknown sex'
    self.test_dir = 'test3'
    self.file_names  = ['test.out', 'test.inf', 'out'  ]
    self.cmd         = 'testrelmatrix data2.par data2.ped >out 2>&1' 
    self.execute()

  def testrelpair(self):
    'relpair test'
    self.test_dir = 'testrel'
    self.file_names  = ['testrelpair.out' ]
    self.cmd         = 'testrelpair dbh.par dbh.dat > testrelpair.out'
    self.execute()

  def testrelsimple(self):
    'relpair_simple test'
    self.test_dir = 'testrelsimple'
    self.file_names  = ['testrelpair_simple.out' ]
    self.cmd         = 'testrelpair_simple mypara myped1and2 > testrelpair_simple.out'
    self.execute()
