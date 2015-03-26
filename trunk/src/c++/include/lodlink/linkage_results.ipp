//============================================================================
// File:      linkage_results.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/25/2 created        -djb
//                                                                          
// Notes:     Inline implementation of linkage results classes.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  lod_ratio_result
//============================================================================
//
inline
lod_ratio_result::lod_ratio_result()
      : perform_test(false)
{}

inline
lod_ratio_result::~lod_ratio_result()
{}

inline double
lod_ratio_result::p_value_u_bound() const
{
  return perform_test ? pow(10, -lod_score()) : QNAN;
}

inline double
lod_ratio_result::chi_sq_stat() const
{
  return perform_test ? hypothesis_result::chi_sq_stat() : QNAN;
}


//============================================================================
// IMPLEMENTATION:  non_ss_lod_ratio_result
//============================================================================
//
inline
non_ss_lod_ratio_result::non_ss_lod_ratio_result()
      : restricted_alt_theta(QNAN), restricted_alt_theta_ub(QNAN),
        relaxed_alt_theta(QNAN), relaxed_alt_theta_ub(QNAN)
{
  var_cov.resize(1, 1, QNAN);
  var_cov_relaxed.resize(1, 1, QNAN);  
}

inline double
non_ss_lod_ratio_result::p_value() const
{
  return  perform_test && ! SAGE::isnan(chi_sq_stat()) ? (1 - ndtr(pow(chi_sq_stat(), .5))) : QNAN;
}

// - Likelihood ratio test is performed if maximum is a local maximum
//   or theta is fixed by maxfun at 0.
//
inline void
non_ss_lod_ratio_result::set_perform_test(const MAXFUN::Results& data)
{
  perform_test = ! (SAGE::isnan(max_checked_value(data)) || bound_at_point_five(data));
}

// - Since the Maxfun documentation does not describe what is meant by
//   fixed 'near' a bound, assume if parameter is fixed at or near a bound
//   and not equal to 0, it is bound to or near .5.  Note: it is assumed direct
//   search is not the last maximization method used, so statuses 9 and 10
//   can be ignored.
//
inline bool
non_ss_lod_ratio_result::bound_at_point_five(const MAXFUN::Results& data)
{
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  const MAXFUN::Parameter&     parameter = parameter_manager.getParameter(AVERAGE);

  MAXFUN::Parameter::ParamTypeEnum  status = parameter.getFinalType();
  double  value = parameter.getFinalEstimate();
  
  return  ((status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND   ||
            status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        ||
            status == MAXFUN::Parameter::IND_FUNC_FIXED_NEAR_BOUND ||
            status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND        ) &&
            value != 0.0);
}

inline void
non_ss_lod_ratio_result::dump(ostream& out) const
{
  out << "\nnon-sex-specific ratio test" << endl;
  hypothesis_result::dump(out);
}

inline void  
non_ss_lod_ratio_result::write_summary(ostream& out) const
{
  //dump(cout);

  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(non_ss_lod_ratio_result::columns.offset()) << ""
  
      << left
      << setw(non_ss_lod_ratio_result::columns.col_w(0))
      << marker
      << setw(non_ss_lod_ratio_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_lod_ratio_result::columns.col_w(1)); 
      write_double(out, restricted_alt_theta); 
  out << setw(non_ss_lod_ratio_result::columns.spc_w(1)) << ""
     
      << setw(non_ss_lod_ratio_result::columns.col_w(2)); 
      write_double(out, relaxed_alt_theta);
  out << setw(non_ss_lod_ratio_result::columns.spc_w(2)) << ""
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(non_ss_lod_ratio_result::columns.col_w(3)); 
      write_double(out, lod_score());
  out << setw(non_ss_lod_ratio_result::columns.spc_w(3)) << ""
  
      << setw(non_ss_lod_ratio_result::columns.col_w(4)); 
      write_double(out, chi_sq_stat());
  out << setw(non_ss_lod_ratio_result::columns.spc_w(4)) << ""
  
      << setprecision(PRC3)
      << setw(non_ss_lod_ratio_result::columns.col_w(5)); 
      write_double(out, p_value(), true);
  out << setw(non_ss_lod_ratio_result::columns.spc_w(5)) << ""
  
      << setw(non_ss_lod_ratio_result::columns.col_w(6)); 
      write_double(out, p_value_u_bound(), true);
  out << setw(non_ss_lod_ratio_result::columns.spc_w(6)) << ""  
      
      << endl;
  
  out.flags(old_flags);
}
 
inline void  
non_ss_lod_ratio_result::write_detail(ostream& out) const
{
 
}

inline void
non_ss_lod_ratio_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}


//============================================================================
// IMPLEMENTATION:  ss_lod_ratio_result
//============================================================================
//
inline
ss_lod_ratio_result::ss_lod_ratio_result()
      : restricted_alt_thetas(QNAN, QNAN), restricted_alt_theta_ubs(QNAN, QNAN),
        relaxed_alt_thetas(QNAN, QNAN), relaxed_alt_theta_ubs(QNAN, QNAN)
{
  var_cov.resize(2, 2, QNAN);
  var_cov_relaxed.resize(2, 2, QNAN);  
}

// - This is not in accordance w. the 3.1 manual, but the way RCE wants it
//   for now.  See SCR 621 for a later refinement.
//
inline double
ss_lod_ratio_result::p_value() const
{
  if(perform_test && ! SAGE::isnan(chi_sq_stat()))
  {
    double  pv = 1;
  
    pv -= .25;
    pv -= .5 * chdtr(1, chi_sq_stat());
    pv -= .25 * chdtr(2, chi_sq_stat());
  
    return  pv;  
  }
  else
  {
    return  QNAN;
  }
}

// - Likelihood ratio test is performed if maximum is a local maximum
//   or both thetas are fixed by maxfun at 0 or one theta is fixed at zero
//   and the partial derivative relative to the other one is zero.
//
inline void
ss_lod_ratio_result::set_perform_test(const MAXFUN::Results& data)
{
  perform_test = ! (SAGE::isnan(max_checked_value(data)) || either_bound_at_point_five(data));
}

// - Since the Maxfun documentation does not describe what is meant by
//   fixed 'near' a bound, assume if parameter is fixed at or near a bound
//   and not equal to 0, it is bound to or near .5.  Note: it is assumed direct
//   search is not the last maximization method used, so statuses 9 and 10
//   can be ignored.
//
inline bool
ss_lod_ratio_result::either_bound_at_point_five(const MAXFUN::Results& data)
{
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  const MAXFUN::Parameter&     male_parameter = parameter_manager.getParameter(MALE);
  const MAXFUN::Parameter&     female_parameter = parameter_manager.getParameter(FEMALE);

  MAXFUN::Parameter::ParamTypeEnum  male_status   = male_parameter.getFinalType();
  MAXFUN::Parameter::ParamTypeEnum  female_status = female_parameter.getFinalType();

  double  male_value   = male_parameter.getFinalEstimate();
  double  female_value = female_parameter.getFinalEstimate();

  bool  male_bound = ((male_status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND   ||
                       male_status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        ||
                       male_status == MAXFUN::Parameter::IND_FUNC_FIXED_NEAR_BOUND ||
                       male_status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND        ) &&
                       male_value != 0.0);

  bool  female_bound = ((female_status == MAXFUN::Parameter::IND_FUNC_FIXED_AT_BOUND   ||
                         female_status == MAXFUN::Parameter::IND_FIXED_AT_BOUND        ||
                         female_status == MAXFUN::Parameter::IND_FUNC_FIXED_NEAR_BOUND ||
                         female_status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND        ) &&
                         female_value != 0.0);

  return  male_bound || female_bound;  
}

inline void
ss_lod_ratio_result::dump(ostream& out) const
{
  out << "\nsex-specific ratio test" << endl;
  hypothesis_result::dump(out);
}

inline void  
ss_lod_ratio_result::write_summary(ostream& out) const
{
  //dump(cout);

    ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(ss_lod_ratio_result::columns.offset()) << ""
  
      << left
      << setw(ss_lod_ratio_result::columns.col_w(0))
      << marker
      << setw(ss_lod_ratio_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_lod_ratio_result::columns.col_w(1)); 
      write_double(out, restricted_alt_thetas.male_theta); 
  out << setw(ss_lod_ratio_result::columns.spc_w(1)) << ""
     
      << setw(ss_lod_ratio_result::columns.col_w(2)); 
      write_double(out, restricted_alt_thetas.female_theta);
  out << setw(ss_lod_ratio_result::columns.spc_w(2)) << ""
  
      << setw(ss_lod_ratio_result::columns.col_w(3)); 
      write_double(out, relaxed_alt_thetas.male_theta); 
  out << setw(ss_lod_ratio_result::columns.spc_w(3)) << ""
     
      << setw(ss_lod_ratio_result::columns.col_w(4)); 
      write_double(out, relaxed_alt_thetas.female_theta);
  out << setw(ss_lod_ratio_result::columns.spc_w(4)) << ""  
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(ss_lod_ratio_result::columns.col_w(5)); 
      write_double(out, lod_score());
  out << setw(ss_lod_ratio_result::columns.spc_w(5)) << ""
  
      << setw(ss_lod_ratio_result::columns.col_w(6)); 
      write_double(out, chi_sq_stat());
  out << setw(ss_lod_ratio_result::columns.spc_w(6)) << ""
  
      << setprecision(PRC3)
      << setw(ss_lod_ratio_result::columns.col_w(7)); 
      write_double(out, p_value(), true);
  out << setw(ss_lod_ratio_result::columns.spc_w(7)) << ""
  
      << setw(ss_lod_ratio_result::columns.col_w(8)); 
      write_double(out, p_value_u_bound(), true);
  out << setw(ss_lod_ratio_result::columns.spc_w(8)) << ""  
      
      << endl;
  
  out.flags(old_flags);
}
 
inline void  
ss_lod_ratio_result::write_detail(ostream& out) const
{

}

inline void
ss_lod_ratio_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}


//============================================================================
// IMPLEMENTATION:  cleves_elston_result
//============================================================================
//
inline
cleves_elston_result::cleves_elston_result()
      : alt_thetas(QNAN, QNAN), alt_theta_ubs(QNAN, QNAN),
        null_thetas(QNAN, QNAN), null_theta_ubs(QNAN, QNAN)
{
  var_cov.resize(2, 2, QNAN);
}

inline double
cleves_elston_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()) ||
     (alt_thetas.male_theta > .5 && alt_thetas.female_theta > .5))
  {
    return  QNAN;
  }
  else
  {
    return  1 - ndtr(pow(chi_sq_stat(), .5));
  }
}

inline void
cleves_elston_result::dump(ostream& out) const
{
  out << "\nCleves-Elston test" << endl;
  hypothesis_result::dump(out);
}

inline void  
cleves_elston_result::write_summary(ostream& out) const
{
  //dump(cout);

  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(cleves_elston_result::columns.offset()) << ""
  
      << left
      << setw(cleves_elston_result::columns.col_w(0))
      << marker
      << setw(cleves_elston_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(cleves_elston_result::columns.col_w(1)); 
      write_double(out, alt_thetas.male_theta); 
  out << setw(cleves_elston_result::columns.spc_w(1)) << ""
     
      << setw(cleves_elston_result::columns.col_w(2)); 
      write_double(out, alt_thetas.female_theta);
  out << setw(cleves_elston_result::columns.spc_w(2)) << ""
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(cleves_elston_result::columns.col_w(3)); 
      write_double(out, lod_score());
  out << setw(cleves_elston_result::columns.spc_w(3)) << ""
  
      << setw(cleves_elston_result::columns.col_w(4)); 
      write_double(out, chi_sq_stat());
  out << setw(cleves_elston_result::columns.spc_w(4)) << ""
  
      << setprecision(PRC3)
      << setw(cleves_elston_result::columns.col_w(5)); 
      write_double(out, p_value(), true);
  out << setw(cleves_elston_result::columns.spc_w(5)) << ""
      
      << endl;
  
  out.flags(old_flags);
}
 
inline void  
cleves_elston_result::write_detail(ostream& out) const
{
  
}

inline void
cleves_elston_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}

