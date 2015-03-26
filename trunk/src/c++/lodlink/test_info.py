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

class parser_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testparser.inf', 'testparser.out']

  def test_parser(self):
    'test parser'
    self.test_dir = 'parser'
    self.cmd         = 'test_parser par ped mld > /dev/null 2>&1'    
    self.execute()


class mle_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testmle.out']

  def test_mle(self):
    'test mle_sub_model'
    self.test_dir = 'mle'
    self.cmd         = 'test_mle > /dev/null 2>&1'    
    self.execute()
    
    
class tcalc_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testtcalc.out']

  def test_tcalc(self):
    'test trans_calculator'
    self.test_dir = 'tcalc'
    self.cmd         = 'test_tcalc > /dev/null 2>&1'    
    self.execute()    


class peeler_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testpeeler.out']

  def test_peeler1(self):
    'simple pedigree, single recombination fraction'
    self.test_dir = 'peeler1'
    self.cmd         = 'test_peeler par ped mld .1 .1 > /dev/null 2>&1'    
    self.execute()    
    
  def test_peeler2(self):
    'simple pedigree, sex_specific recombination fractions'
    self.test_dir = 'peeler2'
    self.cmd         = 'test_peeler par ped mld .2 .3 > /dev/null 2>&1'    
    self.execute()        
    
  def test_peeler3(self):
    'very simple pedigree, sex_specific recombination fractions'
    self.test_dir = 'peeler3'
    self.cmd         = 'test_peeler par ped mld .2 .3 > /dev/null 2>&1'    
    self.execute()        
    
  def test_peeler4(self):
    'equation 8 Cleves/Elston 1997'
    self.test_dir = 'peeler4'
    self.cmd         = 'test_peeler par ped mld .2 .3 > /dev/null 2>&1'
    self.execute()        
    
  def test_peeler5(self):
    'equation 8 Cleves/Elston 1997.  Extended pedigree'
    self.test_dir = 'peeler5'
    self.cmd         = 'test_peeler par ped mld .4 .1 > /dev/null 2>&1'    
    self.execute()
    
    
class mpcalc_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testmpcalc.out']
    
  def test_mpcalc1(self):
    'lodscore should match that of peeler4'
    self.test_dir = 'mpcalc1'
    self.cmd         = 'test_mpcalc par ped mld .2 .3 > /dev/null 2>&1'    
    self.execute()        
    
  def test_mpcalc2(self):
    'same calculation as in peeler2 except alpha = .3'
    self.test_dir = 'mpcalc2'
    self.cmd         = 'test_mpcalc par ped mld .2 .3 .3 > /dev/null 2>&1'    
    self.execute()        
    
  def test_mpcalc_dbh1(self):
    'pedigree 5 of dbh data.  single theta'
    self.test_dir = 'mpcalc_dbh1'
    self.cmd         = 'test_mpcalc par ped mld .3 .3 > /dev/null 2>&1'    
    self.execute()    
    

class max_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['testmax.out']
   
  # - NOTE: LONG RUN TIME
  # 
  def test_max1(self):
    'maximization.  pedigree 5 dbh data.  single theta. LONG RUN TIME.'
    self.test_dir = 'max_dbh1'
    self.cmd      = 'test_max par ped mld false > diagnostic.out 2>&1'    
    self.execute()    
    

class loop_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path      = 'tests'
    self.file_names       = ['screen']
    self.success_expected = 0
    self.cmd              = 'lodlink par ped mld > screen 2>&1'

  def test_loops(self):
    'two loops formed by siblings marrying siblings'
    self.test_dir = 'test_loops'
    self.execute()                    

  def test_marriage_rings(self):
    'one marriage ring'
    self.test_dir = 'test_marriage_rings'
    self.execute()                    


class lodlink_tests (sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['lodlink_analysis1.sum', 'lodlink_analysis2.sum', 
                        'lodlink_analysis1.det', 'lodlink_analysis2.det', 
                        'genome.inf', 'lodlink.inf', 'screen'            ]
    self.cmd         = 'lodlink par ped mld > screen 2>&1'

  def test_lods(self):
    'lod score calculations verified by Cleves-Elston, 1997'
    self.test_dir = 'lods'
    self.execute()

  def test_non_ss_missing(self):
    'missing sex information, non sex specific analysis'
    self.file_names = ['lodlink_analysis1.sum', 'lodlink_analysis1.det',
                       'genome.inf', 'lodlink.inf', 'screen'            ]
    self.test_dir = 'test_non_ss_missing'
    self.execute()                
    
  def test_ss_missing(self):
    'missing sex information, sex specific analysis'
    self.file_names = ['lodlink_analysis1.sum', 'lodlink_analysis1.det',
                       'genome.inf', 'lodlink.inf', 'screen'            ]
    self.test_dir = 'test_ss_missing'
    self.execute()                
    
  def test_x_linkage(self):
    'check test for x_linked markers'
    self.file_names = ['lodlink_analysis1.sum', 'lodlink_analysis1.det',
                       'genome.inf', 'lodlink.inf', 'screen'            ]
    self.test_dir = 'test_x_linkage'
    self.execute()                    

  # - NOTE: LONG RUN TIME.
  #
  def test_sim_ss(self):
    'sex specific tests on simulated data.'
    self.test_dir = 'sim_ss'
    self.delta = .001
    self.execute()           

  # - NOTE: LONG RUN TIME.
  #
  def test_sim_non_ss(self):
    'non sex specific tests on simulated data.'
    self.test_dir = 'sim_non_ss'
    self.delta = .02
    self.execute()                   
 
  def test_sim_ss_mort(self):
    'Mortons test on sex specific simulated data'
    self.file_names = ['lodlink_analysis1.det', 'genome.inf', 'lodlink.inf', 'screen']
    self.test_dir = 'sim_ss_mort'
    self.delta = 0.02
    self.execute()
                    
    
  def test_sim_non_ss_mort(self):
    'Mortons test on non sex specific simulated data'
    self.file_names = ['lodlink_analysis1.det', 'genome.inf', 'lodlink.inf', 'screen']    
    self.test_dir = 'sim_non_ss_mort'
    self.delta = 0.02
    self.execute()


"""
          #'dbh_segreg'      : ('lodlink par ped mld typ', ['lodlink_analysis1.sum', 'lodlink_analysis1.det', 'genome.inf', 'lodlink.inf'],
          #                     'dbh data w. type file from SEGREG'),                                                                                    
          #'sim_segreg'      : ('lodlink par ped mld typ', ['lodlink_analysis1.sum', 'lodlink_analysis1.det', 'genome.inf', 'lodlink.inf'],
          #                     'simulated data w. type file from SEGREG'),                                                                                              
          #'genotypes_non_ss' : ('lodlink par ped mld typ', ['lodlink_analysis1.sum', 'lodlink_analysis1.det', 'genome.inf', 'lodlink.inf'],
          #                     'genotype probabilities using non sex-specific simulated data w. type file from SEGREG'),                                                                                    
          #'genotypes_ss'     : ('lodlink par ped mld typ', ['lodlink_analysis1.sum', 'lodlink_analysis1.det', 'genome.inf', 'lodlink.inf'],
          #                      'genotype probabilities using sex-specific simulated data w. type file from SEGREG')                                                                                                        
"""                                    