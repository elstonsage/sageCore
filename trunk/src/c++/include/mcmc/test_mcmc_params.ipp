// ------------------------------------------------------
// Inline Implementation of test_mcmc_parameters
// ------------------------------------------------------

inline bool
test_mcmc_parameters::is_multipoint() const
{
  return my_multipoint;
}

inline test_mcmc_region_iterator
test_mcmc_parameters::region_begin() const
{
  return my_regions.begin();
}

inline test_mcmc_region_iterator
test_mcmc_parameters::region_end() const
{
  return my_regions.end();
}

inline size_t
test_mcmc_parameters::region_count() const
{
  return my_regions.size();
}

inline void
test_mcmc_parameters::set_multipoint(bool m) 
{
  my_multipoint = m;
}

inline void
test_mcmc_parameters::add_region(const string& s, const string& o)
{
  string t = toUpper(s);

  bool found = false;

  test_mcmc_region_list::iterator i = my_regions.begin();
  for( ; !found && i != region_end(); ++i )
    if( i->name == t )
      found = true;

  if( !found )
    my_regions.push_back(test_mcmc_region_type(t, o));
}

inline void
test_mcmc_parameters::remove_region(const string& s)
{
  string t = toUpper(s);

  test_mcmc_region_list::iterator i = my_regions.begin();
  for( ; i != my_regions.end(); ++i )
  {
    if( i->name == t )
    {
      my_regions.erase(i);
    }
  }
}

