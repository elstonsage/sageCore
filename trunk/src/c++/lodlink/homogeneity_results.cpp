//============================================================================
// File:      homogeneity_results.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   3/4/3 - created.                         djb
//                                                                          
// Notes:     Implementation of homogeneity results classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/homogeneity_results.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  non_ss_mortons_result
//============================================================================
//
header  non_ss_mortons_result::meta;
header  non_ss_mortons_result::columns;
header  non_ss_mortons_result::detail_columns;

void
non_ss_mortons_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_MORTONS_NON_SS_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_LOCUS);
  detail_columns.add_col(_RECOM_LARGE_);
  detail_columns.add_col(_LIKELIHOOD_);
}


//============================================================================
// IMPLEMENTATION:  ss_mortons_result
//============================================================================
//
header  ss_mortons_result::meta;
header  ss_mortons_result::columns;
header  ss_mortons_result::detail_columns;

void
ss_mortons_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_MORTONS_SS_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_LOCUS);
  detail_columns.add_col(_M_RECOM_LARGE_);
  detail_columns.add_col(_F_RECOM_LARGE_);  
  detail_columns.add_col(_LIKELIHOOD_);
}

}
}



