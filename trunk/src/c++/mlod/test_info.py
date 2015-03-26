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


class mlod_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'
    self.file_names = ['out', 'test.det']

  def test_sim1(self):
    'Tests simulated MLOD data'
    self.cmd = 'mlod -p par -d genesis.ped -m trt -l loc -g genome 2>&1 >out'
    self.execute()

  def test_back_message(self):
    "Tests MLOD's production of backwards compatibility error messages"
    self.success_expected = 0
    self.cmd = 'mlod -p par -d ped -m trt -l loc -g genome >out 2>&1 '
    self.file_names = ['out']
    self.execute()

  def test_no_data(self):
    "Tests MLOD's reporting when there is no valid data due to pedigrees being too large"
    self.success_expected = 1
    self.cmd = 'mlod -p Chr01_Rep03.PRM -d Chr01_Rep03 -m rep3_2rec.typ -l rep3freq1.loc -g SageMap.MP >out 2>&1'
    self.file_names = ['out', 'mlod.inf', 'mlod1rep3.det', 'mlod1rep3.sum' ]
    self.execute()

  def test_nan(self):
    "Tests re-scaling of mpoint likelihood, not to report Non-number"
    self.success_expected = 1
    self.cmd = 'mlod -p par -d ped -m trt -l loc -g gen >out 2>&1'
    self.file_names = ['out', 'mlod.inf', 'genome.inf', 'mlod_analysis1.det', 'mlod_analysis1.sum' ]
    self.execute()
