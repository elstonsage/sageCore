//============================================================================
// File:      results.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   12/19/2 - created.                         djb
//                                                                          
// Notes:     Implementation of misc classes.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/results.h"

namespace SAGE
{

namespace LODLINK
{

bool  
converged(const MAXFUN::Results& data)
{
  bool return_value = false;
  
  int  last_error = data.getExitFlag();
  if(0 < last_error && last_error < 4)
  {
    return_value = true;
  }  
  
  return  return_value;
}

// - Determine if a maximum was found by checking last error status for a 
//   convergence.  If a convergence status was not returned,
//   determine if a maximum was found by checking the individual status
//   of each parameter.  If each parameter is fixed to or near a bound or 
//   the partial derivative w. respect to the parameter is near 0, then a
//   maximum was found.
//   Note:  The above implies a Maxfun exit flag is a convergence.
//
//   Return Maxfun_Data.value() if maximum was found.  QNAN, otherwise.
//
double
max_checked_value(const MAXFUN::Results& data)
{
  if(converged(data))
  {
    return  data.getFinalFunctionValue();
  }
  
  bool  max_found = true;
  int  independent_param_count = data.getNumOfVaryingIndependentParams();
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  for(int i = 0; i < independent_param_count; ++i)
  {
    const MAXFUN::Parameter&  parameter = parameter_manager.getParameter(i);
    if(one_non_zero_std_error(data) && parameter.isDerivAvailable() && abs(parameter.getDeriv()) < DERIV_EPS) 
    {
      continue;
    }
    else
    {
      max_found = false;
      break;
    }
  }
  
  return  max_found ? data.getFinalFunctionValue() : QNAN;
}


bool
one_non_zero_std_error(const MAXFUN::Results& data)
{
  bool  one_non_zero = false;
  
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  int  param_count = parameter_manager.getParamCount();
  for(int i = 0; i < param_count; ++i)
  {
    // - Fixed externally to Maxfun.  Per GCW std error will be 0 anyway.
    //
    const MAXFUN::Parameter&  parameter = parameter_manager.getParameter(i);
    if(parameter.getFinalType() == MAXFUN::Parameter::FIXED)
    {
      continue;
    }
    
    if(parameter.isStdErrorAvailable() && parameter.getStdError() != 0)
    {
      one_non_zero = true;
      break;
    } 
  }
  
  return  one_non_zero;
}


//============================================================================
// IMPLEMENTATION:  hypothesis_result
//============================================================================
//
void  
hypothesis_result::set_var_cov(const MAXFUN::Results& data)
{
  const MAXFUN::CovarianceMatrix&  cv_matrix = data.getCovarianceMatrix();
  if(cv_matrix.isAvailable())
  {
    size_t  pc = static_cast<size_t>(cv_matrix.getSize());
    for(size_t row = 0; row < pc; ++row)
    {
      for(size_t col = 0; col < pc; ++col)
      {
        var_cov(row, col) = cv_matrix.getCovariance(row, col);
      }
    }
  }
}

void  
hypothesis_result::write_var_cov(ostream& out) const
{
  if(var_cov.rows() == 0)
  {
    return;           // Nothing to write.
  }

  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ') << left;
  out << setw(_LOCUS.tw()) << marker;
  write_vc_row(out, var_cov, 0);
  out << endl;
      
  for(size_t row = 1; row < var_cov.rows(); ++row)
  {
    out << setw(_LOCUS.tw()) << "";
    write_vc_row(out, var_cov, row);
    out << endl;    
  }
  
  out << endl;
  out.flags(old_flags);  
}
 
void  
hypothesis_result::write_vc_row(ostream& out, const Matrix2D<double>& matrix, size_t row)
{
  assert(row < matrix.rows());

  ios::fmtflags old_flags = out.flags();
  
  out << right << setprecision(PRC4);
  for(size_t col = 0; col < matrix.cols(); ++col)
  {
    out << setw(VAR_COV_SZ);
    write_double(out, matrix(row, col));
  }
  
  out.flags(old_flags);
}


//============================================================================
// IMPLEMENTATION:  non_ss_smiths_result
//============================================================================
//
header  non_ss_smiths_result::meta;
header  non_ss_smiths_result::columns;
header  non_ss_smiths_result::vc_meta;

void
non_ss_smiths_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_SMITHS_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_RECOM_SMALL_);
  columns.add_col(_EST_LINKAGE_PROP_);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  vc_meta.set_offset(_LOCUS.tw() + 2);
  vc_meta.add_col(_SMITHS_NON_SS_VC_);  
}


//============================================================================
// IMPLEMENTATION:  ss_smiths_result
//============================================================================
//
header  ss_smiths_result::meta;
header  ss_smiths_result::columns;
header  ss_smiths_result::vc_meta;

void
ss_smiths_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_SMITHS_);
  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_M_RECOM_SMALL_);
  columns.add_col(_F_RECOM_SMALL_);  
  columns.add_col(_EST_LINKAGE_PROP_);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  vc_meta.set_offset(_LOCUS.tw() + 2);
  vc_meta.add_col(_SMITHS_SS_VC_);    
}


//============================================================================
// IMPLEMENTATION:  non_ss_faraways_result
//============================================================================
//
header  non_ss_faraways_result::meta;
header  non_ss_faraways_result::columns;
header  non_ss_faraways_result::detail_columns;
header  non_ss_faraways_result::vc_meta;

void
non_ss_faraways_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_FARAWAYS_);
  
  // - Summary
  //  
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_RECOM_SMALL_);
  columns.add_col(_EST_LINKAGE_PROP_);
  columns.add_col(_LOD);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  // - Detail
  //
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_TWO_LINE_LOCUS);
  detail_columns.add_col(_POSTERIOR);  
  
  // - Variance-covariance
  //
  vc_meta.set_offset(_LOCUS.tw() + 2);
  vc_meta.add_col(_FARAWAYS_NON_SS_VC_);    
}


//============================================================================
// IMPLEMENTATION:  ss_faraways_result
//============================================================================
//
header  ss_faraways_result::meta;
header  ss_faraways_result::columns;
header  ss_faraways_result::detail_columns;
header  ss_faraways_result::vc_meta;

void
ss_faraways_result::build_headers()
{
  meta.set_offset(_LOCUS.tw());
  meta.add_col(_FARAWAYS_);
  
  // - Summary
  //
  columns.set_underline(U_CHR);
  columns.add_col(_LOCUS);
  columns.add_col(_M_RECOM_SMALL_);
  columns.add_col(_F_RECOM_SMALL_);  
  columns.add_col(_EST_LINKAGE_PROP_);
  columns.add_col(_LOD);
  columns.add_col(_CHI);
  columns.add_col(_P_VALUE);
  
  // - Detail
  //
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_TWO_LINE_LOCUS);
  detail_columns.add_col(_POSTERIOR);
  
  // - Variance-covariance
  //
  vc_meta.set_offset(_LOCUS.tw() + 2);
  vc_meta.add_col(_FARAWAYS_SS_VC_);    
}

}
}
