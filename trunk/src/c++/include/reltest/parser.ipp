// ---------------------------------------------------------------------------
// Inline Implementation of reltest_parser
// ---------------------------------------------------------------------------

inline
const string&
reltest_parser::get_analysis_name() const
{
  return my_analysis_name;
}

inline
const vector<size_t>&
reltest_parser::get_analysis_regions() const
{
  return my_analysis_regions;
}

inline
const vector<bool>&
reltest_parser::get_analysis_pairtypes() const
{
  return my_analysis_pairtypes;
}

inline
const vector<double>&
reltest_parser::get_preset_cutpoints() const
{
  return my_preset_cutpoints;
}

inline
bool
reltest_parser::calculate_cutpoints() const
{
  return my_calculate_cutpoints;
}

inline
bool
reltest_parser::generate_nucfam_output() const
{
  return my_nucfam_out;
}

inline
bool
reltest_parser::generate_detailed_output() const
{
  return my_detailed_out;
}

