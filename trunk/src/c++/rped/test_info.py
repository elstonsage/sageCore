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

class rped_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_rp_info(self):
    'Test the RPED Info classes'
    self.cmd = "test_rp_info >out 2>&1"
    self.file_names = ['out']
    self.execute()

  def test_mpfile1(self):
    'MPfile Test 1'
    self.cmd = "mpfiletest ped.dat '(2X,A2,1X,A2,T18,A1,T8,2A3)' out > mpfiletest.out 2>&1"
    self.file_names = ['mpfiletest.out']
    self.execute()

  def test_mpfile2(self):
    'MPfile Test 2'
    self.cmd = "mpfiletest ped2.dat '(2(A4,1X),T21,A1,T10,2(1X,A4))' out2 > mpfiletest.out 2>&1"
    self.file_names = ['mpfiletest.out']
    self.execute()

  def test_loops(self):
    'MultiPedigree Loop Tests'
    self.relative_path = 'tests/test_loops'
    self.cmd = 'test_loops >/dev/null 2>&1'
    self.file_names = ['loop_test.out']
    self.execute()

  def test_x_linked(self):
    'Test of x-linked marker parsing'
    self.cmd = '../app/test_app params ped >/dev/null 2>&1'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_noped_nopar_nosex(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_noped_nopar_sex(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_noped_par_nosex_1(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_noped_par_nosex_2(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_noped_par_nosex_2(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_noped_par_sex(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_nopar_nosex_1(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_nopar_nosex_2(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_nopar_sex_1(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_nopar_sex_2(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_par_nosex_1(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ped_par_nosex_2(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_ped_par_sex(self):
    self.cmd = '../app/test_app params ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

  def test_par_sex_err(self):
    self.cmd = '../app/test_app params ped >out 2>&1'
    self.file_names = ['ageon.inf', 'out']
    self.success_expected=0
    self.execute()

  def test_category(self):
    self.cmd = '../app/test_app par ped >/dev/null 2>/dev/null'
    self.file_names = ['ageon.inf']
    self.execute()

class marker_list_tests (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test_ml_no_params(self):
    'marker_list without start or end attributes test'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_no_st(self):
    'marker_list without start attribute'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_no_end(self):
    'marker_list without end attribute'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_st_not_in_ped(self):
    'marker_list with start marker not in pedigree file'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_end_not_in_ped(self):
    'marker_list with end marker not in pedigree file'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_end_b4_st(self):
    'marker_list with end marker before the start in the pedigree file'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_p_at_st(self):
    'marker_list with start marker also a standard parameter'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_p_at_end(self):
    'marker_list with end marker also a standard parameter'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_p_in_mid(self):
    'marker_list with middle marker also a standard parameter'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_ml_at_st(self):
    'marker_list with start marker also in another marker list'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_ml_at_end(self):
    'marker_list with end marker also in another marker list'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_overlap_ml_in_mid(self):
    'marker_list with middle marker also in another marker list'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=0
    self.execute()

  def test_ml_valid(self):
    'marker_list valid'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()

  def test_ml_valid_2_ml(self):
    'marker_list valid'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()

  def test_ml_valid_1_marker(self):
    'marker_list valid'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()

  def test_ml_valid_2_ml_1_marker(self):
    'marker_list valid'
    self.cmd = "../app/test_app par ../test_marker_list/ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()

  def test_marker_cov(self):
    'marker_covariate'
    self.cmd = "../app/test_app par ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()

  def test_covariate_list(self):
    'covariate_list'
    self.cmd = "../app/test_app par ped >/dev/null 2>&1"
    self.file_names = ['ageon.inf']
    self.success_expected=1
    self.execute()
