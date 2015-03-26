#ifndef SEGREG_UTILITIES_H
#define SEGREG_UTILITIES_H
//===================================================================
//
//  File:	segreg_utilities.h
//
//  Author:	Stephen Gross
//
//  History:	sag Initial implementation.		Aug 03 2001
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===================================================================


#include <iostream>
#include <iomanip>
#include "app/output_streams.h"
#include "rped/rped.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "mped/mp_utilities.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/sub_model_base.h"
#include "segreg/model.h"
#include "segreg/member_calculator.h"
#include "segreg/polygenic_penetrance_calculator.h"
#include "segreg/SL_calculator.h"
#include "segreg/regressive_peeler.h"
#include "segreg/mlm_peeler.h"
#include "segreg/MlmCorrelationVerifier.h"
#include "segreg/mlm_resid_corr.h"
#include "segreg/mlm_fra_util.h"
#include "segreg/FPMM_peeler.h"
#include "segreg/polygenic_transition_calculator.h"
#include "pedcalc/fra_test_util.h"

namespace SAGE { namespace SEGREG {

/// general testing utility functions

/// The purpose of this class is to provide testing functions for the
/// various calculator classes involved in segreg. The functions are all
/// static, so you don't have to instantiate anything in order to test the
/// calculators.

class segreg_utilities
{
  public:

    typedef list<FPED::SubpedigreeConstPointer>              subped_list;

    /// tests the continuous member calculator.
    static int test_continuous_MC(const FPED::Multipedigree & RMP, const model & mod, bool asc);

    /// tests the binary member calculator.
    static int test_binary_MC(const FPED::Multipedigree & RMP, const model & mod, bool asc);

    /// tests the onset member calculator.
    static int test_onset_MC(const FPED::Multipedigree & RMP, const model & mod, bool ascer);

    static int test_binary_fam_resid
                  (const FPED::Multipedigree& RMP,
                   const model&                  mod,
                   bool                          ascer);
                   
    /// test_PENETRANCE(...) tests the (regressive) penetrance calculator.
    static int test_PENETRANCE(const FPED::Multipedigree & RMP,
      const model & mod, bool ascer);

    /// tests the polygenic_penetrance_calculator
    static int test_polygenic_penetrance(const FPED::Multipedigree & RMP,
      const model & mod, bool ascer);

    /// test_POLYGENIC_TRANSITION
    static int test_POLYGENIC_TRANSITION(const FPED::Multipedigree & RMP,
      const model & mod);

    /** test_A(...) tests the SL calculator for Regressive models. */
    static int test_Regressive(const FPED::Multipedigree & RMP, const model & mod);
    static int test_FPMM(const FPED::Multipedigree & RMP, const model & mod);
    static int test_ASCER_CALC(const FPED::Multipedigree & RMP, const model & mod);
    static bool almost_equals(double x, double y);

    static int test_SEGREG_CALCULATOR(const FPED::Multipedigree & RMP, const model & mod);
    static int test_MLM              (const FPED::Multipedigree & RMP, const model & mod, bool ascer);
    static void   analysis_controller  (const FPED::Multipedigree & RMP, APP::Output_Streams& os);
    static string fp                   (double d, size_t w, size_t p);

    static void test_prev_calculation(const model& mod);
                    
    /// test the mlm_pair_correlations pair generation algorithm
    static void test_mlm_pairs(const PedigreeDataSet& RMP, const model& mod);
    
    /// test the mlm_resid_correlations pair generation algorithm
    static void test_mlm_correlations(const PedigreeDataSet& RMP, const model& mod);
};

// End namespace
}}

#endif
