
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

class unr_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'hapapp_analysis1.det', 'hapapp_analysis1.sum', 'screen']

  def test_unr0(self):
    """
    Purpose:  test population frequency estimation.
              test generation of all diplotypes. 
              test calculation of diplotype probabilities
    
    Basis:    hand calculation.
              inspection. 
    """
    self.test_dir = 'test_unr0'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()
    
  def test_unr1(self):
    """
    Purpose:  test population frequency estimation.
    
    Basis:    Python prototype, multi/estimate, which was tested against Arlequin. 
    """
    self.test_dir = 'test_unr1'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()    
    
  def test_unr_region(self):
    """
    Purpose:  test region specification.
              test generation of all possible diplotypes.
              
    Basis:    inspection. 
    """
    self.test_dir = 'test_unr_region'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()    
    
  def test_unr_sub_pops(self):
    """
    Purpose:  test assignment of data members to subpopulations.
    
    Basis:    inspection. 
    """
    self.test_dir = 'test_unr_sub_pops'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()            
    
  def test_unr_x_linkage(self):
    """
    Purpose:  test program w. x-linked loci (unrelated individuals).
    
    Basis:    Hand calculation and inspection. 
    """
    self.test_dir = 'test_unr_x_linkage'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                
    
  
class fam_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'hapapp_analysis1.det', 'hapapp_analysis1.sum', 'screen']
  
  def test_fam_ex2(self):
    """
    Note:     Katrina's example 2.
    
    Purpose:  test use of family information to determine possible diplotypes.
              test assignment of family representatives.
    
    Basis:    hand calculation and inspection. 
    """
    self.test_dir = 'test_fam_ex2'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                
    
  def test_fam_merlin2(self):
    """
    Purpose:  test use of family information to determine possible diplotypes.
    
    Basis:    Merlin output and inspection. 
    """
    self.test_dir = 'test_fam_merlin2'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                    
    
  def test_fam_merlin3(self):
    """
    Purpose:  test use of family information to determine possible diplotypes.
    
    Basis:    Merlin output and inspection. 
    """
    self.test_dir = 'test_fam_merlin3'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()   
    
  def test_fam_merlin4(self):
    """
    Purpose:  test use of family information to determine possible diplotypes.
    
    Basis:    Merlin output and inspection. 
    """
    self.test_dir = 'test_fam_merlin4'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                            

  def test_fam_remap1(self):
    """
    Purpose:  test program response to an invalid analysis block.
    
    Basis:    Inspection. 
    """
    self.test_dir = 'test_fam_remap1'
    self.file_names  = ['decipher.inf', 'screen']    
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                                         
    
  def test_fam_remap2(self):
    """
    Purpose:  test modification of genotype elimination to eliminate remapping
              when subpedigree has all missing data at a marker.
    
    Basis:    Inspection. 
    """
    self.test_dir = 'test_fam_remap2'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                                             
    
  def test_fam_remap3(self):
    """
    Purpose:  test modification of genotype elimination to eliminate remapping
              when subpedigree has all missing data at a marker.
    
    Basis:    Inspection. 
    """
    self.test_dir = 'test_fam_remap3'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                                                 

  def test_partitioning(self):
    """
    Purpose:  Test selection of family reps and subpopulation members with two
              partition variables.

    Basis:    Inspection.  Diplotype combinations not verified.
    """
    self.test_dir = 'test_partitioning'
    self.file_names += ['hapapp_analysis2.det', 'hapapp_analysis2.sum',
                        'hapapp_analysis3.det', 'hapapp_analysis3.sum' ]
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()
    
class sim_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'decipher_analysis1.det', 'decipher_analysis1.sum', 'screen']
  
  def test_fam_sim(self):
    """
    Purpose:  test use of family information in determining haplotype
              frequency estimates.
    
    Basis:    simulation. 
    """
    self.test_dir = 'test_fam_sim'
    self.cmd      = 'decipher -p par -d ped > screen 2>&1'
    self.execute()
    
    
class fpool_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'decipher_analysis1.sum', 'screen']
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'

  def test_inconsistency(self):
    """
    Purpose:  test handling of Mendelian inconsistencies. 
    
    Basis:    inspection. 
    """
    self.test_dir = 'test_inconsistency'
    self.file_names = ['decipher.inf', 'rep.sum', 'rep.det', 'fpool.sum', 'fpool.det', 'screen']
    self.execute()                            
  
  def test_fpool_sim1(self):
    """
    Purpose:  test haplotype frequency estimates using founder pool analysis unit 
              w simulated data.
    
    Basis:    simulation, Merlin, Haplore. 
    """
    self.test_dir = 'test_fpool_sim1'
    self.execute()                        
    
  def test_mixed_data(self):
    """
    Purpose:  force founder pool option into singleton and family_rep mode for part of  
              the data.
    
    Basis:    inspection, test_fpool_1B. 
    """
    self.test_dir = 'test_mixed_data'
    self.file_names = ['decipher.inf', 'mixed.det', 'screen']
    self.cmd = 'decipher -p mixed_par -d mixed_ped -g gen > screen 2>&1'
    self.execute()                            
    
    
class fpool_combs_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher_analysis1.det', 'screen']    
    self.cmd         = 'decipher -p par -d ped -g gen > screen 2>&1'    
  
  def test_fpool_1B(self):
    """
    Purpose:  test determination of founder pool haplotype combinations. 
    
    Basis:    examples submitted by Rob Igo (see .ppt files in .../decipher/tests). 
    """
    self.test_dir = 'test_fpool_1B'
    self.execute()                    
    
  def test_fpool_1J(self):
    """
    Purpose:  test determination of founder pool haplotype combinations. 
    
    Basis:    examples submitted by Rob Igo (see .ppt files in .../decipher/tests). 
    """
    self.test_dir = 'test_fpool_1J'
    self.execute()                        
    
  def test_fpool_2B(self):
    """
    Purpose:  test determination of founder pool haplotype combinations. 
    
    Basis:    examples submitted by Rob Igo (see .ppt files in .../decipher/tests).
              2nd combination in 2nd row should be 0-0-0 / 0-1-1 as confirmed by
              Rob, however.
 
    """
    self.test_dir = 'test_fpool_2B'
    self.execute()                

  def test_fpool_incons(self):
    """
    Purpose:  test determination of founder pool haplotype combinations when
              one locus contains a Mendelian inconsistancy. 
    
    Basis:    inspection.
 
    """
    self.test_dir = 'test_fpool_incons'
    self.execute()                



class pool_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'analysis1.det', 'analysis1.sum', 'screen']
    self.cmd      = 'decipher -p par -d ped > screen 2>&1'    
  
  def test_pool_vs_unr0(self):
    """
    Purpose:  test use of pool information in determining all possible
              haplotype combinations, most likely haplotype combinations and 
              haplotype frequency estimates.
    
    Basis:    test_unr0.  This test uses the same data, but in pool format. 
    """
    self.test_dir = 'test_pool_vs_unr0'
    self.execute()                    
   

class pool_comb_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['analysis1.det', 'screen']    
    self.cmd         = 'decipher -p par -d ped  > screen 2>&1'    
  
  def test_pool1(self):
    """
    Purpose:  test determination of pool haplotype combinations. 
    
    Basis:    Inspection. 
    """
    self.test_dir = 'test_pool1'
    self.execute()                             
    
  def test_pool2(self):
    """
    Purpose:  test determination of pool haplotype combinations. 
    
    Basis:    Inspection. 
    """
    self.test_dir = 'test_pool2'
    self.execute()                                 
    
    
class pool_sim_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'analysis1.det', 'analysis1.sum', 'screen']
    self.cmd      = 'decipher -p par -d ped > screen 2>&1'    
  
  def test_pool_sim1(self):
    """
    Purpose:  test use of pool information in determining haplotype
              frequency estimates.
    
    Basis:    simulation (see .../decipher/tests/pool_sim_script). 
    """
    self.test_dir = 'test_pool_sim1'
    self.execute()                    
    
  def test_pool_sim2(self):
    """
    Purpose:  test use of pool information in determining haplotype
              frequency estimates.
    
    Basis:    simulation (see .../decipher/tests/pool_sim_script).
              Comparison to the program, EHP.R. 
    """
    self.test_dir = 'test_pool_sim2'
    self.file_names  = ['decipher.inf', 'analysis1.sum', 'screen']    
    self.execute()


class marker_filtration_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'decipher_analysis1.sum', 'screen']
  
  def test_maf1(self):
    """
    Purpose:  test marker filtration by minor allele frequency.  This test
              also illustrates the priority of the marker locus description
              file over the parameter file pedigree block in determining marker
              order. 
    
    Basis:    Inspection.  Haplotype frequency estimates not verified. 
    """
    self.test_dir = 'test_maf1'
    self.cmd      = 'decipher -p par -d ped -l mld > screen 2>&1'
    self.execute()                        

    
class block_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'decipher_analysis1.sum', 'screen']
  
  def test_four_gamete1(self):
    """
    Purpose:  test four gamete rule implementation.
    
    Basis:    Haploview (see README file in test directory). 
    """
    self.test_dir = 'test_four_gamete1'
    self.cmd      = 'decipher -p par -d ped > screen 2>&1'
    self.execute()    
    
  def test_sliding_window1(self):
    """
    Purpose:  test sliding window implementation.
    
    Basis:    inspection (NOTE: haplotype frequencies not verified). 
    """
    self.test_dir = 'test_sliding_window1'
    self.cmd      = 'decipher -p par -d ped > screen 2>&1'
    self.execute()        
    
  def test_ld2(self):
    """
    Purpose:  test ld block determination.  Also provides a test of the dump file option.
    
    Basis:    LD calculation verified vs haploview.  Block determination verified by inspection
              (NOTE: haplotype frequencies not verified). 
    """
    self.file_names  = ['decipher.inf', 'decipher_analysis1.sum', 'decipher_analysis1.det', 'decipher_analysis1.dmp', 'screen']    
    self.test_dir    = 'test_ld2'
    self.delta       = .01
    self.cmd         = 'decipher -p par -d ped -l mld > screen 2>&1'
    self.execute()            
    

""" These tests take too long to run on a routine basis    
class lr_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'
    self.file_names  = ['decipher.inf', 'decipher_analysis1.det', 'decipher_analysis1.sum', 'screen']
  
  def test_lr3(self):
    
    #Purpose:  test likelihood ratio test and calculation of empirical p-values.
    
    #Basis:    simulation and asymptotic vs empirical p-values.

    #Note:     should not be part of routine testing as it takes 2 hours to run.
    #          6-14-6.  Modified so that it now takes 4 1/2 minutes to run.
 
    
    self.test_dir = 'test_lr3'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute() 

  def test_lr4(self):
    
    #Purpose:  test likelihood ratio test and calculation of empirical p-values.
    
    #Basis:    simulation and asymptotic vs empirical p-values. 
    
    #Note:     should not be part of routine testing as it takes 20 hours to run!
    #          'exp' files not up to date with regard to output format changes.
    
    
    self.test_dir = 'test_lr4'
    self.cmd      = 'decipher -p par -d ped -g gen > screen 2>&1'
    self.execute()                        
"""