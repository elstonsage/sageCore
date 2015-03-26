////////////////////////////////////////////////////////////////////////////
//             Implementation of ARP_base_analysis(Inline)                //
////////////////////////////////////////////////////////////////////////////

inline
void
ARP_base_analysis::set_parameters(const lodpal_parameters& p)
{
  invalidate();
  my_parameters = p;
}

inline
void
ARP_base_analysis::set_previous_parameters(const lodpal_parameters& pp)
{
  my_previous_parameters = pp;
}

inline
void
ARP_base_analysis::set_good_bound(bool b)
{
  my_good_bound = b;
}

inline
bool
ARP_base_analysis::is_good_bound() const
{
  return my_good_bound;
}

inline
lodpal_parameters&
ARP_base_analysis::parameters()
{
  return my_parameters;
}

inline
const lodpal_parameters&
ARP_base_analysis::parameters() const
{
  return my_parameters;
}

inline
lodpal_parameters&
ARP_base_analysis::previous_parameters()
{
  return my_previous_parameters;
}

inline
const lodpal_parameters&
ARP_base_analysis::previous_parameters() const
{
  return my_previous_parameters;
}

inline
lodpal_pairs&
ARP_base_analysis::pairs_info()
{
  return my_pairs;
}

inline
const lodpal_pairs&
ARP_base_analysis::pairs_info() const
{
  return my_pairs;
}

inline
const RelativePairs&
ARP_base_analysis::relative_pairs() const
{
  return my_pairs.relative_pairs();
}

inline
size_t
ARP_base_analysis::valid_parameter_count() const
{
   return my_valid_parameter_count;
}

inline
bool
ARP_base_analysis::built() const
{
   if(!parameters().valid())
     return false;

   return my_built;
}

inline
double
ARP_base_analysis::get_lod_mm() const
{
  return my_ld_mm;
}

inline
double
ARP_base_analysis::get_lod_mf() const
{
  return my_ld_mf;
}

inline
double
ARP_base_analysis::get_lod_ff() const
{
  return my_ld_ff;
}

inline
lodpal_result&
ARP_base_analysis::get_lodpal_result()
{
  return my_result;
}

inline
const lodpal_result&
ARP_base_analysis::get_lodpal_result() const
{
  return my_result;
}

inline
const marker_type& 
ARP_base_analysis::current_marker() const
{
  return my_current_marker;
}

inline
void
ARP_base_analysis::set_current_marker(const marker_type& m)
{
  my_current_marker = m;
}
