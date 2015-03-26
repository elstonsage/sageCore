#ifndef LODLINK_HOMOGENEITY_TESTS_H
#define LODLINK_HOMOGENEITY_TESTS_H
//============================================================================
// File:      homogeneity_tests.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/13/3 created         -djb
//                                                                          
// Notes:     Defines classes for calculating and writing results of 
//            linkage homogeneity tests.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "maxfun/sub_model.h"
#include "lodlink/instructions.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"
#include "lodlink/homogeneity_results.h"
#include "lodlink/tasks.h"

using std::vector;
using std::string;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    non_ss_mortons_test
//                                                                          
//  Purpose:  represent Morton's test for linkage homogeneity w. a single
//            recombination fraction.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class non_ss_mortons_test : public task
{
    friend void do_task_calculations<non_ss_mortons_test, non_ss_mortons_result>(non_ss_mortons_test&);

  public:
    non_ss_mortons_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);

    void  announce_start() const;
    void  calculate();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, non_ss_mortons_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, non_ss_mortons_result& result);  
    
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
};

//----------------------------------------------------------------------------
//  Class:    ss_mortons_test
//                                                                          
//  Purpose:  represent Morton's test for linkage homogeneity w. separate
//            male and female recombination fractions.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class ss_mortons_test : public task
{
  friend void do_task_calculations<ss_mortons_test, ss_mortons_result>(ss_mortons_test&);

  public:
    ss_mortons_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);

    void  announce_start() const;
    void  calculate();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, ss_mortons_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, ss_mortons_result& result);  
    
    void  write_summary(ostream& out) const;
    void  write_detail(ostream& out) const;
};

#include "lodlink/homogeneity_tests.ipp"

}
}

#endif
