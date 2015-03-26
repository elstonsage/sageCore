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

class genibd_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['exact_multi.CHR3.ibd', 'exact_multi.X.ibd',
                        'single.CHR3.ibd', 'single.X.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.cmd         = 'genibd -p par -d ped -l loc -g gen >out 2>&1'

  def test_sib1(self):
    'simple nuclear family test : only full sib pairs'
    self.test_dir = 'test_sib1'
    self.execute()
    
  def test_sib2(self):
    'simple extended pedigree test : full sib & grand parental pairs'
    self.test_dir = 'test_sib2'
    self.execute()

  def test_rel1(self):
    'all 5 different relative pair types'
    self.test_dir = 'test_rel1'
    self.execute()

  def test_all_pair(self):
    'All pair option test'
    self.file_names  = ['exact_multi.CHR3.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.execute()

  def test_ex_and_sim(self):
    'Tests intervals with both exact and simulation methods, which previously had some data missing in output'
    self.file_names  = ['multi.CHR3.ibd']
    self.execute()

  def test_rel1_sim(self):
    'all 5 different relative pair types, using simulation'
    self.test_dir = 'test_rel1_sim'
    self.file_names = ['sim_single.CHR3.ibd', 'sim_multi.CHR3.ibd']
    self.execute()

  def test_rel2(self):
    'multipoint with only one marker : should give the same result as singlepoint from test_rel1'
    self.test_dir = 'test_rel2'
    self.file_names  = ['Analysis_1.CHR3_1.ibd', 'Analysis_4.X1.ibd',
                        'Analysis_2.CHR3_2.ibd', 'Analysis_5.X2.ibd',
                        'Analysis_3.CHR3_3.ibd', 'Analysis_6.X3.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.execute()

  def test_split(self):
    'pedigree split test'
    self.test_dir = 'test_split'
    self.file_names  = ['exact_multi.CHR3.ibd', 'exact_multi.CHR4.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.execute()

  def test_split2(self):
    'pedigree is splitted at the filtering stage without split_pedigrees=true.'
    self.test_dir = 'test_split2'
    self.file_names  = ['multi.CH_17.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.execute()

  # test needed when sex-inference algo is done.
  def test_nosex(self):
    'The sex of a parent or both parents is missing. -> infer or assign arbitrary sex'
    self.test_dir = 'test_nosex'
    self.file_names  = ['hyp10_3.CHR_1.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
                        
    self.delta = 0.01
    self.execute()
    
  def test_spoint(self):
    'Tests singlepoint calculation on some GAW pedigrees'
    self.test_dir = 'test_spoint'
    self.file_names = ['out', 'single.CH22.ibd' ]
    self.execute()

  def test_sim_sp(self):
    'Tests simulation singlepoint calculation on some GAW pedigrees'
    self.test_dir = 'test_sim_sp'
    self.file_names = ['single.CH22.ibd' ]
    self.execute()

  def test_x(self):
    'Tests basic x-linked case which might exibit double-Y'
    self.file_names = ['single.CH22.ibd' ]
    self.file_names  = ['exact_multi.X.ibd',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.cmd         = 'genibd par ped loc gen >out 2>&1'
    self.execute()

  def test_ibd_state(self):
    'IBD state output'
    self.test_dir = 'test_ibd_state'
    self.file_names  = ['exact_multi.CHR3.ibd', 'exact_multi.CHR3.state',
                        'out', 'genibd.inf', 'genome.inf'  ]
    self.execute()
    
  # take too long to run, so use only when needed.
  #def test_snp(self):
  #  'GAW14 affy-SNP test'
  #  self.test_dir = 'test_snp'
  #  self.file_names  = ['affy.CHROMOSOME1.ibd',
  #                      'out', 'genibd.inf', 'genome.inf'  ]
  #  self.execute()
    
