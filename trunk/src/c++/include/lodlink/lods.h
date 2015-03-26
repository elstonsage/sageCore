#ifndef LODLINK_LODS_H
#define LODLINK_LODS_H
//============================================================================
// File:      lods.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/13/3 created         -djb
//                                                                          
// Notes:     Defines classes for calculating lod scores.
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
#include "lodlink/instructions.h"
#include "lodlink/tasks.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"
#include "lodlink/lods_results.h"

using std::vector;
using std::string;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    non_ss_lods
//                                                                          
//  Purpose:  calculates lod scores for specified non-sex specific 
//            recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
class non_ss_lods : public task
{
    friend void do_task_calculations<non_ss_lods, non_ss_lods_result>(non_ss_lods&);

  public:
    non_ss_lods(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    
    void  announce_start() const;
    void  calculate();

  private:  
    void  calculate_alt(size_t trait_index, size_t marker_index, non_ss_lods_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, non_ss_lods_result& result);
    
    void  write_summary(ostream& out) const;
      void  write_summary_headers(ostream& out) const;
      void  write_columns_header(ostream& out) const;
    void  write_detail(ostream& out) const;
      static void  write_detail_header(ostream& out);
};

//----------------------------------------------------------------------------
//  Class:    ss_lods_test
//                                                                          
//  Purpose:  calculates lod scores for specified sex specific
//            recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
class ss_lods : public task
{
    friend void do_task_calculations<ss_lods, ss_lods_result>(ss_lods&);

  public:
    ss_lods(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    
    void  announce_start() const;
    void  calculate();

  private:
    void  calculate_alt(size_t trait_index, size_t marker_index, ss_lods_result& result);
    void  calculate_null(size_t trait_index, size_t marker_index, ss_lods_result& result);
    
    void  write_summary(ostream& out) const;
      void  write_summary_headers(ostream& out) const;
      void  write_columns_header(ostream& out) const;
    void  write_detail(ostream& out) const;
      static void  write_detail_header(ostream& out);
};

#include "lodlink/lods.ipp"

}
}

#endif
