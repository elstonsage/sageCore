// ======================================
// Inline functions of sim_ibd_analysis
// ======================================
  
inline
cerrorstream
sim_ibd_analysis::get_errors() const
{
  return errors;
}

inline
cerrorstream
sim_ibd_analysis::set_errors(SAGE::cerrorstream& s)
{
  errors = s;
  return errors;
}

inline bool
sim_ibd_analysis::built() const
{
  return my_built;
}

inline bool
sim_ibd_analysis::valid() const
{
  return my_valid;
}

inline SAGE::IBD*
sim_ibd_analysis::ibd_adaptor() const
{
  return my_ibds;
}

