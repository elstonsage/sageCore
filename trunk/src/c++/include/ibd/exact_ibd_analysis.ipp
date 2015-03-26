//---------------------------------------------------------------------------
// Inline Implementation of exact_ibd_analysis
//---------------------------------------------------------------------------

inline
cerrorstream
exact_ibd_analysis::get_errors() const
{
  return errors;
}

inline
cerrorstream
exact_ibd_analysis::set_errors(cerrorstream& s)
{
  errors = s;
  my_ldata.set_errors(s);

  return errors;
}

inline
bool
exact_ibd_analysis::built() const
{
  return my_built;
}

inline
bool
exact_ibd_analysis::valid() const
{
  return my_valid;
}

inline
IBD*
exact_ibd_analysis::ibd_adaptor() const
{
  return my_ibds;
}

//inline
//IBD*
//exact_ibd_analysis::sib_ibd_adaptor() const
//{
//  return my_sib_ibds;
//}
