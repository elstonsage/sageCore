//============================================================================
// File:      homogeneity_results.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/13/3 created        -djb
//                                                                          
// Notes:     Inline implementation of homogeneity results classes.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  non_ss_mortons_result
//============================================================================
//
inline
non_ss_mortons_result::non_ss_mortons_result()
      : group_theta_ub(QNAN), null_theta(QNAN), null_theta_ub(QNAN)
{}

inline double
non_ss_mortons_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return  chdtrc(group_results.size() - 1, chi_sq_stat());
  }
}

inline void  
non_ss_mortons_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(non_ss_mortons_result::columns.offset()) << ""
  
      << left
      << setw(non_ss_mortons_result::columns.col_w(0))
      << marker
      << setw(non_ss_mortons_result::columns.spc_w(0)) << ""
      
      << right
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(non_ss_mortons_result::columns.col_w(1)); 
      write_double(out, chi_sq_stat());
  out << setw(non_ss_mortons_result::columns.spc_w(1)) << ""
  
      << setprecision(PRC3)
      << setw(non_ss_mortons_result::columns.col_w(2)); 
      write_double(out, p_value(), true);
  out << setw(non_ss_mortons_result::columns.spc_w(2)) << ""
      << endl;
  
  out.flags(old_flags);
}

inline void  
non_ss_mortons_result::write_detail(ostream& out) const
{

}

inline void
non_ss_mortons_result::write_group_detail(ostream& out, const string& group) const
{
  std::map<string, pair<double, double> >::const_iterator  g_iter = group_results.find(group);
  assert(g_iter != group_results.end());

  ios::fmtflags old_flags = out.flags();

  out << setfill(' ');
  out << setw(non_ss_mortons_result::detail_columns.offset()) << ""
  
      << left
      << setw(non_ss_mortons_result::detail_columns.col_w(0))
      << marker
      << setw(non_ss_mortons_result::detail_columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_mortons_result::detail_columns.col_w(1)); 
  write_double(out, g_iter->second.first);
  out << setw(non_ss_mortons_result::detail_columns.spc_w(1)) << ""  
  
      << resetiosflags(ios::floatfield)  << setprecision(PRC3)
      << setw(non_ss_mortons_result::detail_columns.col_w(2)); 
  write_double(out, g_iter->second.second);      
  
  out << endl;
  
  out.flags(old_flags);  
}

inline void
non_ss_mortons_result::write_vc_matrix(ostream& out) const
{

}


//============================================================================
// IMPLEMENTATION:  ss_mortons_result
//============================================================================
//
inline
ss_mortons_result::ss_mortons_result()
      : group_theta_ubs(QNAN, QNAN), null_thetas(QNAN, QNAN), 
        null_theta_ubs(QNAN, QNAN)
{}

inline double
ss_mortons_result::p_value() const
{
  if(SAGE::isnan(chi_sq_stat()))
  {
    return  QNAN;
  }
  else
  {
    return  chdtrc(2 * (group_results.size() - 1), chi_sq_stat());
  }
}

inline void  
ss_mortons_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << setw(ss_mortons_result::columns.offset()) << ""
  
      << left
      << setw(ss_mortons_result::columns.col_w(0))
      << marker
      << setw(ss_mortons_result::columns.spc_w(0)) << ""
      
      << right
      << resetiosflags(ios::floatfield)  << setprecision(PRC2)
      << setw(ss_mortons_result::columns.col_w(1)); 
      write_double(out, chi_sq_stat());
  out << setw(ss_mortons_result::columns.spc_w(1)) << ""
  
      << setprecision(PRC3)
      << setw(ss_mortons_result::columns.col_w(2)); 
      write_double(out, p_value(), true);
  out << setw(ss_mortons_result::columns.spc_w(2)) << ""
      
      << endl;
  
  out.flags(old_flags);
}

inline void  
ss_mortons_result::write_detail(ostream& out) const
{

}

inline void
ss_mortons_result::write_group_detail(ostream& out, const string& group) const
{
  std::map<string, pair<theta_pair, double> >::const_iterator  g_iter = group_results.find(group);
  assert(g_iter != group_results.end());

  ios::fmtflags old_flags = out.flags();

  out << setfill(' ');
  out << setw(ss_mortons_result::detail_columns.offset()) << ""
  
      << left
      << setw(ss_mortons_result::detail_columns.col_w(0))
      << marker
      << setw(ss_mortons_result::detail_columns.spc_w(0)) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_mortons_result::detail_columns.col_w(1)); 
  write_double(out, g_iter->second.first.male_theta);
  out << setw(ss_mortons_result::detail_columns.spc_w(1)) << ""  
  
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_mortons_result::detail_columns.col_w(2)); 
  write_double(out, g_iter->second.first.female_theta);
  out << setw(ss_mortons_result::detail_columns.spc_w(2)) << ""    
  
      << resetiosflags(ios::floatfield)  << setprecision(PRC3)
      << setw(ss_mortons_result::detail_columns.col_w(3)); 
  write_double(out, g_iter->second.second);      
  
  out << endl;
  
  out.flags(old_flags);  
}

inline void
ss_mortons_result::write_vc_matrix(ostream& out) const
{

}


