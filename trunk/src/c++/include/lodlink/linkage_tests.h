#ifndef LODLINK_LINKAGE_TESTS_H
#define LODLINK_LINKAGE_TESTS_H
//============================================================================
// File:      linkage_tests.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/25/2 created         -djb
//                                                                          
// Notes:     Defines classes for calculating and writing results of 
//            linkage tests.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <vector>
#include <cmath>
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "maxfun/sub_model.h"
#include "fped/fped.h"
#include "lodlink/instructions.h"
#include "lodlink/tasks.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"
#include "lodlink/output.h"
#include "lodlink/linkage_results.h"

using std::vector;
using std::string;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    non_ss_lod_ratio_test
//                                                                          
//  Purpose:  represents a lod ratio test for linkage w. a non-sex specific
//            recombination fraction.
//                                                                          
//----------------------------------------------------------------------------
//
class non_ss_lod_ratio_test : public task
{
  public:

    friend void do_task_calculations<non_ss_lod_ratio_test, non_ss_lod_ratio_result>(non_ss_lod_ratio_test&);

    non_ss_lod_ratio_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    ~non_ss_lod_ratio_test();
    
    void  announce_start() const;
    void  calculate();

  private:  
    void  calculate_alt(size_t trait_index, size_t marker_index, non_ss_lod_ratio_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, non_ss_lod_ratio_result& result);
    void  calculate_relaxed_theta(size_t trait_index, size_t marker_index, 
                                  non_ss_lod_ratio_result& result, bool point_five);
    static bool  est_at_point_five(const MAXFUN::Results& data);
    
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
};

//----------------------------------------------------------------------------
//  Class:    ss_lod_ratio_test
//                                                                          
//  Purpose:  represents a lod ratio test for linkage w. sex specific
//            recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
class ss_lod_ratio_test : public task
{
    friend void do_task_calculations<ss_lod_ratio_test, ss_lod_ratio_result>(ss_lod_ratio_test&);

  public:
    ss_lod_ratio_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    ~ss_lod_ratio_test();
    
    void  announce_start() const;
    void  calculate();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, ss_lod_ratio_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, ss_lod_ratio_result& result);
    void  calculate_relaxed_theta(size_t trait_index, size_t marker_index, 
                                  ss_lod_ratio_result& result, bool male_point_five, bool female_point_five);
    void  calculate_relaxed_theta_bound(size_t trait_index, size_t marker_index, 
                                  ss_lod_ratio_result& result);
    void  calculate_relaxed_theta_compl(size_t trait_index, size_t marker_index, 
                                  ss_lod_ratio_result& result);
    static bool  male_est_at_point_five(const MAXFUN::Results& data);
    static bool  female_est_at_point_five(const MAXFUN::Results& data);
    
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
};


//----------------------------------------------------------------------------
//  Class:    cleves_elston_test
//                                                                          
//  Purpose:  represents cleves_elston test for linkage.
//                                                                          
//----------------------------------------------------------------------------
//
class cleves_elston_test : public task
{
    friend void do_task_calculations<cleves_elston_test, cleves_elston_result>(cleves_elston_test&);

  public:
    cleves_elston_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    ~cleves_elston_test();
    
    void  announce_start() const;
    void  calculate();
    
  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, cleves_elston_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, cleves_elston_result& result);
    
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
};


#include "lodlink/linkage_tests.ipp"

}
}

#endif
