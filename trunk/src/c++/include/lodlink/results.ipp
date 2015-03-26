//============================================================================
// File:      results.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

// - Calculate lod score from two doubles representing the natural
//   logs of the alternative and null likelihoods respectively.
//
inline double
lod_score(double alt_ln_like, double null_ln_like)
{
  return  (alt_ln_like - null_ln_like) / log(10.0);
}


//============================================================================
// IMPLEMENTATION:  task_result
//============================================================================
//
inline
task_result::task_result()
{}

inline
task_result::~task_result()
{}


//============================================================================
// IMPLEMENTATION:  hypothesis_result
//============================================================================
//
inline
hypothesis_result::hypothesis_result()
      : null_ln_like(QNAN), alt_ln_like(QNAN)
{}

inline
hypothesis_result::~hypothesis_result()
{}

// - Twice the difference of the natural logs of the likelihoods.
//
inline double
hypothesis_result::chi_sq_stat() const
{
  double  ls = lod_score();
  
  if(ls >= 0)
  {
    return 2 * log(10.0) * ls;
  }
  else
  {
    return QNAN;
  }
}

inline double
hypothesis_result::lod_score() const
{
  return LODLINK::lod_score(alt_ln_like, null_ln_like);
}

inline void
hypothesis_result::dump(ostream& out) const
{
  out << "marker  " << marker
      << "  alt ln likelihood  " << alt_ln_like 
      << "  null ln likelihood  " << null_ln_like << endl;
  
}

//============================================================================
// IMPLEMENTATION:  subpedigree_posterior
//============================================================================
//
inline
subpedigree_posterior::subpedigree_posterior(const string& mn, double l)
      : member_name(mn), posterior(l)
{}

inline bool
subpedigree_posterior::operator==(const subpedigree_posterior& right) const
{
  return  member_name == right.member_name;
}

//============================================================================
// IMPLEMENTATION:  pedigree_posterior
//============================================================================
//
inline
pedigree_posterior::pedigree_posterior(const string& pn, double l)
      : pedigree_name(pn), posterior(l)
{}

inline bool
pedigree_posterior::operator==(const pedigree_posterior& right) const
{
  return  pedigree_name == right.pedigree_name;
}


//============================================================================
// IMPLEMENTATION:  non_ss_alt_result
//============================================================================
//
inline
non_ss_alt_result::non_ss_alt_result()
      : alt_ln_like(QNAN), alt_theta(QNAN), alt_theta_ub(QNAN), 
        alpha(QNAN)
{}


//============================================================================
// IMPLEMENTATION:  non_smiths_faraways_result
//============================================================================
//
inline
non_ss_smiths_faraways_result::non_ss_smiths_faraways_result()
      : alt_theta(QNAN), alt_theta_ub(QNAN), 
        null_theta(QNAN), null_theta_ub(QNAN), alpha(QNAN)
{
  var_cov.resize(2, 2, QNAN);
}


//============================================================================
// IMPLEMENTATION:  ss_alt_result
//============================================================================
//
inline
ss_alt_result::ss_alt_result()
      : alt_ln_like(QNAN), alt_thetas(QNAN, QNAN), alt_theta_ubs(QNAN, QNAN),
        alpha(QNAN)
{}


//============================================================================
// IMPLEMENTATION:  ss_smiths_faraways_result
//============================================================================
//
inline
ss_smiths_faraways_result::ss_smiths_faraways_result()
      : alt_thetas(QNAN, QNAN), alt_theta_ubs(QNAN, QNAN),
        null_thetas(QNAN, QNAN), null_theta_ubs(QNAN, QNAN), alpha(QNAN)
{
  var_cov.resize(3, 3, QNAN);
}



//============================================================================
// IMPLEMENTATION:  non_ss_faraways_result
//============================================================================
//
inline double
non_ss_faraways_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return   .5 * (1 - pow(chdtr(1, chi_sq_stat()), 2));
  }
}

inline double  
non_ss_faraways_result::get_subpedigree_posterior(const string& pedigree, const string& member) const
{
  vector<pedigree_posterior>::const_iterator  ped_iter;
  ped_iter = find(posteriors.begin(), posteriors.end(), pedigree_posterior(pedigree));
  
  if(ped_iter != posteriors.end())
  {
    vector<subpedigree_posterior>::const_iterator  subped_iter;  
    subped_iter = find(ped_iter->sub_posteriors.begin(),
                       ped_iter->sub_posteriors.end(),
                       subpedigree_posterior(member));
                       
    if(subped_iter != ped_iter->sub_posteriors.end())
    {
      return  subped_iter->posterior;
    }
  }
  
  return  QNAN;
}

inline void  
non_ss_faraways_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(non_ss_faraways_result::columns.offset()) << ""
  
      << left
      << setw(non_ss_faraways_result::columns.col_w(0))
      << marker
      << setw(non_ss_faraways_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_faraways_result::columns.col_w(1)); 
      write_double(out, alt_theta); 
  out << setw(non_ss_faraways_result::columns.spc_w(1)) << ""
     
      << setw(non_ss_faraways_result::columns.col_w(2)); 
      write_double(out, alpha);
  out << setw(non_ss_faraways_result::columns.spc_w(2)) << ""
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(non_ss_faraways_result::columns.col_w(3)); 
      write_double(out, lod_score());
  out << setw(non_ss_faraways_result::columns.spc_w(3)) << ""
  
      << setw(non_ss_faraways_result::columns.col_w(4)); 
      write_double(out, chi_sq_stat());
  out << setw(non_ss_faraways_result::columns.spc_w(4)) << ""
  
      << setprecision(PRC3)
      << setw(non_ss_faraways_result::columns.col_w(5)); 
      write_double(out, p_value(), true);
  out << setw(non_ss_faraways_result::columns.spc_w(5)) << ""
  
      << endl;
  
  out.flags(old_flags);
}
  
inline void  
non_ss_faraways_result::write_detail(ostream& out) const
{

}

inline void  
non_ss_faraways_result::write_family_detail(ostream& out, const string& pedigree, 
                                                          const string& member) const
{
  ios::fmtflags old_flags = out.flags();

  out << setfill(' ');
  out << setw(non_ss_faraways_result::detail_columns.offset()) << ""
  
      << left
      << setw(non_ss_faraways_result::detail_columns.col_w(0))
      << marker
      << setw(non_ss_faraways_result::detail_columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_faraways_result::detail_columns.col_w(1)); 
      
  double  posterior = get_subpedigree_posterior(pedigree, member);
  
  if(! (SAGE::isnan(posterior) || SAGE::isnan(alpha)))
  {
    // - Ott 1991, eq. 9.12, p206
    //
    // - Addition to alpha to eliminate insignificant differences
    //   due to numerical limitations.
    //
    if(posterior > alpha + 100 * numeric_limits<double>::epsilon())
    {
      write_double(out, posterior);
    }
    else
    {
      write_double(out, QNAN);
    }
    
    // #define POST_DEBUG
    #ifdef POST_DEBUG
    out << scientific << setprecision(6)
        << "   " << "difference  " << posterior - alpha 
        << "   " << "required difference   " << 100 * numeric_limits<double>::epsilon(); 
    #endif    
  }
  else
  {
    write_double(out, QNAN);
  }

  out << endl;
  
  out.flags(old_flags);
}

inline void
non_ss_faraways_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}


//============================================================================
// IMPLEMENTATION:  ss_faraways_result
//============================================================================
//
inline double
ss_faraways_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return   1 - pow(chdtr(2, chi_sq_stat()), 2);
  }
}

inline double  
ss_faraways_result::get_subpedigree_posterior(const string& pedigree, const string& member) const
{
  vector<pedigree_posterior>::const_iterator  ped_iter;
  ped_iter = find(posteriors.begin(), posteriors.end(), pedigree_posterior(pedigree));
  
  if(ped_iter != posteriors.end())
  {
    vector<subpedigree_posterior>::const_iterator  subped_iter;  
    subped_iter = find(ped_iter->sub_posteriors.begin(),
                       ped_iter->sub_posteriors.end(),
                       subpedigree_posterior(member));
                       
    if(subped_iter != ped_iter->sub_posteriors.end())
    {
      return  subped_iter->posterior;
    }
  }
  
  return  QNAN;
}

inline void  
ss_faraways_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(ss_faraways_result::columns.offset()) << ""
  
      << left
      << setw(ss_faraways_result::columns.col_w(0))
      << marker
      << setw(ss_faraways_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_faraways_result::columns.col_w(1)); 
      write_double(out, alt_thetas.male_theta); 
  out << setw(ss_faraways_result::columns.spc_w(1)) << ""
  
      << setw(ss_faraways_result::columns.col_w(2)); 
      write_double(out, alt_thetas.female_theta); 
  out << setw(ss_faraways_result::columns.spc_w(2)) << ""  
     
      << setw(ss_faraways_result::columns.col_w(3)); 
      write_double(out, alpha);
  out << setw(ss_faraways_result::columns.spc_w(3)) << ""
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(ss_faraways_result::columns.col_w(4)); 
      write_double(out, lod_score());
  out << setw(ss_faraways_result::columns.spc_w(4)) << ""
  
      << setw(ss_faraways_result::columns.col_w(5)); 
      write_double(out, chi_sq_stat());
  out << setw(ss_faraways_result::columns.spc_w(5)) << ""
  
      << setprecision(PRC3)
      << setw(ss_faraways_result::columns.col_w(6)); 
      write_double(out, p_value(), true);
  out << setw(ss_faraways_result::columns.spc_w(6)) << ""
  
      << endl;
  
  out.flags(old_flags);
}
  
inline void  
ss_faraways_result::write_detail(ostream& out) const
{

}

inline void  
ss_faraways_result::write_family_detail(ostream& out, const string& pedigree, 
                                                      const string& member) const
{
  ios::fmtflags old_flags = out.flags();

  out << setfill(' ');
  out << setw(ss_faraways_result::detail_columns.offset()) << ""
  
      << left
      << setw(ss_faraways_result::detail_columns.col_w(0))
      << marker
      << setw(ss_faraways_result::detail_columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_faraways_result::detail_columns.col_w(1));
      
  double  posterior = get_subpedigree_posterior(pedigree, member);
  
  if(! (SAGE::isnan(posterior) || SAGE::isnan(alpha)))
  {
    // - Ott 1991, eq. 9.12, p206
    //
    // - Addition to alpha to eliminate insignificant differences
    //   due to numerical limitations.
    //    
    if(posterior > alpha + 100 * numeric_limits<double>::epsilon())
    {
      write_double(out, posterior);
    }
    else
    {
      write_double(out, QNAN);
    }
    
    // #define POST_DEBUG
    #ifdef POST_DEBUG
    out << scientific << setprecision(6)
        << "   " << "difference  " << posterior - alpha 
        << "   " << "required difference   " << 100 * numeric_limits<double>::epsilon(); 
    #endif
  }
  else
  {
    write_double(out, QNAN);
  }  
  
  out << endl;
  
  out.flags(old_flags);
}

inline void
ss_faraways_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}



//============================================================================
// IMPLEMENTATION:  non_ss_smiths_result
//============================================================================
//
inline double
non_ss_smiths_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return   .5 * chdtrc(1, chi_sq_stat());
  }
}

inline void  
non_ss_smiths_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(non_ss_smiths_result::columns.offset()) << ""
  
      << left
      << setw(non_ss_smiths_result::columns.col_w(0))
      << marker
      << setw(non_ss_smiths_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_smiths_result::columns.col_w(1)); 
      write_double(out, alt_theta); 
  out << setw(non_ss_smiths_result::columns.spc_w(1)) << ""
     
      << setw(non_ss_smiths_result::columns.col_w(2)); 
      write_double(out, alpha);
  out << setw(non_ss_smiths_result::columns.spc_w(2)) << ""
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(non_ss_smiths_result::columns.col_w(3)); 
      write_double(out, chi_sq_stat());
  out << setw(non_ss_smiths_result::columns.spc_w(3)) << ""
  
      << setprecision(PRC3)
      << setw(non_ss_smiths_result::columns.col_w(4)); 
      write_double(out, p_value(), true);
  out << setw(non_ss_smiths_result::columns.spc_w(4)) << ""
  
      << endl;
  
  out.flags(old_flags);
}
  
inline void  
non_ss_smiths_result::write_detail(ostream& out) const
{

}

inline void
non_ss_smiths_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}


//============================================================================
// IMPLEMENTATION:  ss_smiths_result
//============================================================================
//
inline double
ss_smiths_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return   chdtrc(2, chi_sq_stat());
  }
}

inline void  
ss_smiths_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(ss_smiths_result::columns.offset()) << ""
  
      << left
      << setw(ss_smiths_result::columns.col_w(0))
      << marker
      << setw(ss_smiths_result::columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_smiths_result::columns.col_w(1)); 
      write_double(out, alt_thetas.male_theta); 
  out << setw(ss_smiths_result::columns.spc_w(1)) << ""
  
      << setw(ss_smiths_result::columns.col_w(2)); 
      write_double(out, alt_thetas.female_theta); 
  out << setw(ss_smiths_result::columns.spc_w(2)) << ""  
     
      << setw(ss_smiths_result::columns.col_w(3)); 
      write_double(out, alpha);
  out << setw(ss_smiths_result::columns.spc_w(3)) << ""
      
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(ss_smiths_result::columns.col_w(4)); 
      write_double(out, chi_sq_stat());
  out << setw(ss_smiths_result::columns.spc_w(4)) << ""
  
      << setprecision(PRC3)
      << setw(ss_smiths_result::columns.col_w(5)); 
      write_double(out, p_value(), true);
  out << setw(ss_smiths_result::columns.spc_w(5)) << ""
  
      << endl;
  
  out.flags(old_flags);
}
  
inline void  
ss_smiths_result::write_detail(ostream& out) const
{

}

inline void
ss_smiths_result::write_vc_matrix(ostream& out) const
{
  write_var_cov(out);
}

