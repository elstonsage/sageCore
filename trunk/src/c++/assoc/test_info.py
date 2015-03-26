#============================================================================
#
#	THIS COPY SPECIALIZED FOR ASSOC by Stephen Gross 10 Apr 2003
#
# File:      test_info.py
#                                                                          
# Author:    Dan Baechle
#                                                                          
# History:   Created 7/31/2                                             
#                                                                          
# Notes:     For use by test scripts (see src/c++/test_scripts).
#
#            While in parent of the test directory, use runall to run all tests.
#
#            While in parent of the test directory, use runtest <test directory> 
#            to run a specific test.
#
#            While in test directory, use newexps to make the current
#            test output files the expected files.
#                                                                          
# Copyright (c) 2002 R.C. Elston                                           
# All Rights Reserved                                                    
#============================================================================

import sagetest

class fast_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  # - Added 6-14-7. djb
  #
  def test_sample_desc(self):
    """
    purpose:  test sample description that appears in output files.
    
    basis:    hand calculation.  Only sample description and 
              standardization of residuals verified.
    
    Results for model foo_SEX_CODE changed when program was changed to use null
    estimates as starting values for test model.  P-values close to one in
    both cases, however.
    """
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    self.file_names = [ 'result.sum', 'result_Baseline_null.res' ]
    self.execute()

  """
  def test_binary(self):
    'Binary test'
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    self.delta=0.001
    #self.epsilon=0.007
    self.file_names = [ 'binaryoutput.sum', 'binaryoutput.det', 'out' ]
    self.execute()
  """

  def test_cov_stand(self):
    'Covariate standardizationout'
    self.cmd = 'assoc -p par -d ped >out 2>&1'
    #self.epsilon=0.0001
    self.delta=0.001
    self.file_names = [ 'result.det', 'result.sum', 'out' ]
    self.execute()

  def test_batch(self):
    """
    Coefficient estimates are different for foo_SEX_CODE test case on galton,
    sage_macg5 and sage_cygwin than on marker.  2-16-8.
    
    Model foo_sbpd_rand10 test case has a huge value for the derivative of the sex_code
    coefficient on beastgate.  Also final likelihood for the test case of this model
    is much different on beasgate than on marker.  2-16-8
    """    
    'Batch mode test'
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    #self.epsilon=0.0001
    self.file_names = [ 'assoc.det', 'assoc.sum', 'out' ]
    self.execute()

  def test_trans(self):
    """
    On Beastgate there are differences in somes p-values, se's and derivatives
    between debug and release builds.  1-4-8.    
    """
    'Testing transformation'
    self.cmd = 'assoc -p par -d ped >out 2>&1'
    #self.epsilon=0.015
    self.delta=0.02
    self.file_names = [ 'assoc.inf', 'logsbp.det', 'logsbp.sum', 'out' ]
    self.execute()
    
  def test_nonpolygenic_models(self):
    'Assoc tests'
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    self.epsilon = 0.000001
    self.delta   = 0.005
    self.file_names =  \
                          ['none.det', 'none.sum',
                           'me.det', 'me.sum',
                           'fe.det', 'fe.sum',
                           'feme.det', 'feme.sum',
                           'se.det', 'se.sum',
                           'seme.det', 'seme.sum',
                           'sefe.det', 'sefe.sum',
                           'sefeme.det', 'sefeme.sum', 'out' ]

    self.execute()    


class slow_tests(sagetest.SAGE_Test):

  def setUp(self):
    self.common_path = 'tests'

  def test_polygenic_models(self):
    'Assoc polygenic tests'
    # Note: expected results changed significantly for pe____fe__.det when
    # step-parent offspring correlation no longer calculated (trac ticket #1779).
    # Per RCE, This is because stopping criteria changed w. elimination of this
    # dependent parameter, and the new results should be the correct ones.  Removal
    # or addition of convergence warnings also occured with some models on one
    # or more platforms.
    #
    # Note: expected results changed significantly for pe____fe__.det when
    # correlation logic was corrected (trac ticket #1797).  Removal
    # or addition of convergence warnings also occured with some models on one
    # or more platforms.
    
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    #self.epsilon = 0.001
    self.delta   = 0.01
    self.file_names =  \
                          ['pe.det', 'pe.sum',
                           'peme.det', 'peme.sum',
                           'pefe.det', 'pefe.sum',
                           'pefeme.det', 'pefeme.sum',
                           'pese.det', 'pese.sum',
                           'peseme.det', 'peseme.sum',
                           'pesefe.det', 'pesefe.sum',
                           'pesefeme.det', 'pesefeme.sum', 'out' ]

    self.execute()

  """
  def test_effect_removal(self):
    
    Purpose:  Test removal of variance effects fixed to a bound by MAXFUN.
    
    Basis:    Simulated data supplied by Courtney.
    
    self.cmd = 'assoc -p par -d ped >out 2>&1'
    self.delta = .01
    self.file_names = [ 'assoc.inf', 'p2.det', 'p2.sum', 'out' ]
    self.execute()
    """
    
  """
  def test_user_effects(self):
    
    purpose:  test user-specified variance components.
    
    basis:    tested against model using family effect (see family_rep5.det) which 
              specified exactly the same clusters.  Test data is Courtney's simulation,
              FP_N_NUC_RAND_0005_sample.txt
    
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    self.file_names = [ 'assoc.det', 'assoc.inf', 'out' ]
    self.execute()
    """
    
  def test_user_effects2(self):
    """
    purpose:  test user-specified variance components.
    
    basis:    the user effect in this case should exactly mimic
              the built in family effect.  Estimations of the mean (intercept)
              and the total variance match hand calculations of these quantities
              for the data.
    """
    self.cmd = 'assoc -p par -d ped > out 2>&1'
    self.file_names = [ 'user_fe.det', 'fe.det' ]
    self.execute()

