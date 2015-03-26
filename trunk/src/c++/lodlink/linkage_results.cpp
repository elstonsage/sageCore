//============================================================================
// File:      linkage_results.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   3/3/3 - created.                         djb
//                                                                          
// Notes:     Implementation of linkage results classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/linkage_results.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  lod_ratio_result
//============================================================================
//
void  
lod_ratio_result::set_var_cov_relaxed(const MAXFUN::Results& data)
{
  const MAXFUN::CovarianceMatrix&  vc_matrix = data.getCovarianceMatrix();

  size_t  pc = static_cast<size_t>(vc_matrix.getSize());
  for(size_t row = 0; row < pc; ++row)
  {
    for(size_t col = 0; col < pc; ++col)
    {
      var_cov_relaxed(row, col) = vc_matrix.getCovariance(row, col);
    }
  }
}

void  
lod_ratio_result::write_var_cov(ostream& out) const
{
  if(var_cov.rows() == 0)
  {
    return;           // Nothing to write.
  }

  assert(var_cov.rows() == var_cov_relaxed.rows());

  const size_t  matrix_sw = 20;

  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ') << left;
  out << setw(_LOCUS.tw()) << marker;
  write_vc_row(out, var_cov, 0);
  out << setw(matrix_sw) << "";
  write_vc_row(out, var_cov_relaxed, 0);
  out << endl;
      
  for(size_t row = 1; row < var_cov.rows(); ++row)
  {
    out << setw(_LOCUS.tw()) << "";
    write_vc_row(out, var_cov, row);
    out << setw(matrix_sw) << "";
    write_vc_row(out, var_cov_relaxed, row);    
    out << endl;    
  }
  
  out << endl;
  out.flags(old_flags);  
}


//============================================================================
// IMPLEMENTATION:  non_ss_lod_ratio_result
//============================================================================
//
header  non_ss_lod_ratio_result::meta;
header  non_ss_lod_ratio_result::columns;
header  non_ss_lod_ratio_result::vc_meta;

void
non_ss_lod_ratio_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_MLE_AVE_);
  meta.add_col(_LOD_SCORE_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_RECOM_SMALL_);
  columns.add_col(_RECOM_LARGE_);
  columns.add_col(_LOD);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  columns.add_col(_P_VALUE_UB_);
  
  vc_meta.set_offset(_LOCUS.tw());
  vc_meta.add_col(_NON_SS_LOD_SCORE_VC_);
}


//============================================================================
// IMPLEMENTATION:  ss_lod_ratio_result
//============================================================================
//
header  ss_lod_ratio_result::meta;
header  ss_lod_ratio_result::columns;
header  ss_lod_ratio_result::vc_meta;

void
ss_lod_ratio_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_MLE_SEX_SPEC_);
  meta.add_col(_LOD_SCORE_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_M_RECOM_SMALL_);
  columns.add_col(_F_RECOM_SMALL_);
  columns.add_col(_M_RECOM_LARGE_);
  columns.add_col(_F_RECOM_LARGE_);
  columns.add_col(_LOD);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  columns.add_col(_P_VALUE_UB_);
  
  vc_meta.set_offset(_LOCUS.tw());
  vc_meta.add_col(_SS_LOD_SCORE_VC_);  
}


//============================================================================
// IMPLEMENTATION:  cleves_elston_result
//============================================================================
//
header  cleves_elston_result::meta;
header  cleves_elston_result::columns;
header  cleves_elston_result::vc_meta;

void
cleves_elston_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_CLEVES_ELSTON_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_M_RECOM_LARGE_);
  columns.add_col(_F_RECOM_LARGE_);
  columns.add_col(_LOD);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  vc_meta.set_offset(_LOCUS.tw() + 2);
  vc_meta.add_col(_CLEVES_ELSTON_VC_);
}

}
}



