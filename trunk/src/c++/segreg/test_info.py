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

class parser_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testparser.inf', 'testparser.out']
    self.cmd         = 'test_parser par ped 2>/dev/null >/dev/null '

  def test01(self):
    'misc parameters and mean, transformation and frequency sub-models'
    self.cmd = self.cmd + '1'
    self.execute()

  def test02(self):
    'covariate, resid, variance and fpmm sub-models'
    self.cmd = self.cmd + '2'
    self.execute()

  def test03(self):
    'various sub-models'
    self.cmd = self.cmd + '3'
    self.execute()

  def test04(self):
    'sub-model dump() functions'
    self.cmd = self.cmd + '4'
    self.execute()

  def test06(self):
    'specification of transmission analysis'
    self.cmd = self.cmd + '0'
    self.execute()

  def test07(self):
    'onset sub-model'
    self.cmd = self.cmd + '5'
    self.execute()

  def test08(self):
    'ascertainment sub-model'
    self.cmd = self.cmd + '6'
    self.execute()

  def test09(self):
    'set_two_to_dom() and set_two_to_rec() functions in mean and variance sub-models'
    self.cmd = self.cmd + '7'
    self.execute()

  def test11(self):
    'prevalence estimate sub-model'
    self.cmd = self.cmd + '9'
    self.execute()

  def test12(self):
    'susceptibility sub-model'
    self.cmd = self.cmd + '10'
    self.execute()

  def test13(self):
    'primary trait type'
    self.cmd = self.cmd + '11'
    self.execute()

  def test14(self):
    'covariate sub-models'
    self.cmd = self.cmd + '12'
    self.execute()

  def test15(self):
    'covariates from prevalence sub-models'
    self.cmd = self.cmd + '13'
    self.execute()

  def test16(self):
    'age of onset model options'
    self.cmd = self.cmd + '0'
    self.execute()

class type_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'
    self.file_names = ['out']
  
  def test_type_description(self):
    'Test of the basic type description mechanism'
    self.cmd = 'test_type_description >out 2>&1'
    self.execute()

class numeric_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'
    self.file_names = ['out']

  def test_continuous(self):
    'Tests of the continuous (regressive) MC and penetrance objects'
    self.cmd = 'test_segreg -p par -d fam >out'
    self.execute()

  def test_binary(self):
    'MLM tests (binary member calculator, mlm_peeler'
    self.cmd = 'test_segreg -p binary.par -d fam >out'
    self.execute()

  def test_onset(self):
    'onset member calculator tests'
    self.cmd = 'test_segreg -p onset.par -d fam >out'
    self.execute()

  def test_analysis_sex_effect(self):
    'Testing the analysis sex effect functions'
    self.cmd = "test_segreg -p par -d ped >out 2>&1"
    self.execute()

  def test_mlm_resid(self):
    'Testing of MLM residual correlations computations'
    self.delta       = 0.01
    self.epsilon     = 0.001
    self.cmd         = 'segreg -p par -d fam >out'
    self.file_names  = ["segreg_analysis1.det", "segreg_analysis2.det" ]
    self.execute()

  def test_binary_out(self):
    'MLM tests .typ option'
    self.cmd = 'segreg -p binary.par -d fam >out'
    self.file_names = ["segreg_analysis1.typ"]
    self.execute()

  def test_onset_out(self):
    'Testing age onset model 0 polyenic loci'
    self.cmd = 'segreg -p par -d fam >out'
    self.file_names = ["segreg_analysis1.det"]
    self.execute()

class system_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_iter(self):
    'iterative model tests'
    self.cmd         = 'segreg -p par -d ped >out'
    self.file_names  = ['out', "t1.2seg.det", "t1.com.det", "t1.3general.det",
                        "t1.3seg.det", "t2.2seg.det", "t2.com.det",
                        "t2.3general.det", "t2.3seg.det"]
    self.delta       = 0.05
    self.epsilon     = 0.001
    self.execute()

  def test_mito(self):
    'Testing of mitohondrial DNA'
    self.delta       = 0.01
    self.epsilon     = 0.001
    self.cmd         = 'segreg -p par -d ped >out'
    self.file_names  = ["one_mean_t2.det", "two_mean_t2.det", "three_mean_t2.det" ]
    self.execute()

  def test_transm_no_ped_struct(self):
    'Testing of transmission models missing pedigree structure'
    self.delta       = 0.01
    self.epsilon     = 0.001
    self.cmd         = 'segreg -p par -d fam >out 2>&1'
    self.file_names  = ["out", "segreg.inf" ]
    self.execute()

  def test_three_add(self):
    'Testing of three_add model'
    self.delta       = 0.02
    self.epsilon     = 0.001
    self.cmd         = 'segreg -p par -d ped >out'
    self.file_names  = ["t1.3addmendel.det", "out" ]
    self.execute()

