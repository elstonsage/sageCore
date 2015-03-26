// ------------------------------------------------------
// Inline Implementation of genibd_parameters
// ------------------------------------------------------

inline const string&
genibd_parameters::title() const
{
  return my_title;
}

inline const string&
genibd_parameters::output() const
{
  return my_output.size() ? my_output : my_title;
}

inline bool
genibd_parameters::is_multipoint() const
{
  return my_multipoint;
}

inline bool
genibd_parameters::scan_interval() const
{
  return my_scan_interval;
}

inline double
genibd_parameters::interval_distance() const
{
  return my_interval_distance;
}

inline bool
genibd_parameters::allow_loops() const
{
  return my_loops;
}

inline bool
genibd_parameters::output_ibd_state() const
{
  return my_output_ibd_state;
}

inline choice_type
genibd_parameters::allow_simulation() const
{
  return my_simulation;
}

inline choice_type
genibd_parameters::allow_family_splitting() const
{
  return my_family_split;
}

inline size_type
genibd_parameters::max_exact_size() const
{
  return my_exact_size;
}

inline pair_category_type
genibd_parameters::pair_category() const
{
  return my_pair_type;
}

inline genibd_region_iterator
genibd_parameters::region_begin() const
{
  return my_regions.begin();
}

inline genibd_region_iterator
genibd_parameters::region_end() const
{
  return my_regions.end();
}

inline size_type
genibd_parameters::region_count() const
{
  return my_regions.size();
}

inline
mcmc_parameters&
genibd_parameters::get_sim_parameters()
{
  return my_sim_parameters;
}
             
inline
const mcmc_parameters&
genibd_parameters::get_sim_parameters() const
{
  return my_sim_parameters;
}

inline void
genibd_parameters::set_title(const string& s) 
{
  my_title = s;
}

inline void
genibd_parameters::set_output(const string& s) 
{
  my_output = s;
}

inline void
genibd_parameters::set_multipoint(bool m) 
{
  my_multipoint = m;
}

inline void
genibd_parameters::set_scan_interval(bool b) 
{
  my_scan_interval = b;
}

inline void
genibd_parameters::set_interval_distance(double d) 
{
  my_interval_distance = d;
}

inline void
genibd_parameters::set_allow_loops(bool b) 
{
  my_loops = b;
}

inline void
genibd_parameters::set_output_ibd_state(bool b) 
{
  my_output_ibd_state = b;
}

inline void
genibd_parameters::set_allow_simulation(choice_type b) 
{
  my_simulation = b;
}

inline void
genibd_parameters::set_allow_family_splitting(choice_type b) 
{
  my_family_split = b;
}

inline void
genibd_parameters::set_max_exact_size(long l) 
{
  my_exact_size = l;
}

inline void
genibd_parameters::set_pair_category(pair_category_type t) 
{
  my_pair_type = t;
}

inline void
genibd_parameters::add_region(const string& s, const string& o)
{
  string t = toUpper(s);

  bool found = false;

  genibd_region_list::iterator i = my_regions.begin();
  for( ; !found && i != region_end(); ++i )
    if( i->name == t )
      found = true;

  if( !found )
    my_regions.push_back(genibd_region_type(t, o));
}

inline void
genibd_parameters::remove_region(const string& s)
{
  string t = toUpper(s);

  genibd_region_list::iterator i = my_regions.begin();
  for( ; i != my_regions.end(); ++i )
  {
    if( i->name == t )
    {
      my_regions.erase(i);
    }
  }
}
