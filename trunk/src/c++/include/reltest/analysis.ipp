//---------------------------------------------------------------------------
// Inline Implementation of reltest_analysis
//---------------------------------------------------------------------------

inline
const reltest_parser*
reltest_analysis::get_parser() const
{
  return my_parser;
}

inline
const RPED::RefMultiPedigree*
reltest_analysis::get_multipedigree() const
{
  return my_multipedigree;
}

inline
const vector<putative_pair>&
reltest_analysis::get_ptt_pairs(putative_type p) const
{
  return my_ptt_pairs[p];
}

inline
const vector<putative_pair>&
reltest_analysis::get_misc_pairs(putative_type p) const
{
  return my_misc_pairs[p];
}

inline
const putative_type&
reltest_analysis::get_current_pairtype() const
{
  return my_current_pairtype;
}

inline
double
reltest_analysis::get_chrom_length(size_t c) const
{
  return my_chrom_length[c];
}

inline
double
reltest_analysis::get_total_genome_length() const
{
  return my_total_genome_length;
}

inline
double
reltest_analysis::get_total_map_points() const
{
  return my_total_map_points;
}

inline
double
reltest_analysis::get_cutpoints(cutpoint_type c) const
{
  return my_cutpoints[c];
}

inline
double
reltest_analysis::get_adjusted_cutpoints(cutpoint_type c) const
{
  return my_adjusted_cutpoints[c];
}

inline
double
reltest_analysis::get_AMIC() const
{
  return my_AMIC;
}

inline
double
reltest_analysis::get_Var_Yj() const
{
  return my_Var_Yj;
}

inline
double
reltest_analysis::get_Var_Yjp() const
{
  return my_Var_Yjp;
}

inline
const vector<solution_type>&
reltest_analysis::get_solution_Yj() const
{
  return my_solution_Yj;
}

inline
const vector<solution_type>&
reltest_analysis::get_solution_Yjp() const
{
  return my_solution_Yjp;
}

inline
size_t
reltest_analysis::get_picked_Yj() const
{
  return my_picked_Yj;
}

inline
size_t
reltest_analysis::get_picked_Yjp() const
{
  return my_picked_Yjp;
}

inline
double
reltest_analysis::get_mu_Yj() const
{
  return my_mu_Yj;
}

inline
double
reltest_analysis::get_mu_Yjp() const
{
  return my_mu_Yjp;
}

inline
double
reltest_analysis::get_Yj_mean() const
{
  return my_Yj_mean;
}

inline
double
reltest_analysis::get_standard_error() const
{
  return my_standard_error;
}
