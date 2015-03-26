inline bool
IBD::has_pedigree(const std::string &ped)
{
  error_t e;
  get_pair(ped, "", "", e);
  if(e == bad_ind_both)
    return true;
  else
    return false;
}

inline bool
IBD::use_pair(const std::string &ped, const std::string &i1,
              const std::string &i2, IBD::error_t &e) const
{
  id_pair p = get_pair(ped,i1,i2,e);
  if( e != no_error )
    return use_pair( p.first, p.second );
  else
    return false;
}

inline bool
IBD::valid_pair(const mem_pointer i1, const mem_pointer i2) const
{
  return valid_pair( pair_index(i1, i2) );
}

inline bool
IBD::invalidate_pair(const mem_pointer i1, const mem_pointer i2) const
{
  return invalidate_pair( pair_index(i1, i2) );
}

inline size_t
IBD::add_pair(const std::string &ped, const std::string &i1,
              const std::string &i2, error_t &e)
{
  id_pair p = get_pair(ped, i1, i2, e);
  if( e == no_error )
    return add_pair(p.first, p.second);
  return (size_t) -1;
}

inline size_t
IBD::pair_index(const std::string &ped, const std::string &i1,
                const std::string &i2, error_t &e) const
{
  id_pair p = get_pair(ped, i1, i2, e);
  if( e == no_error )
    return pair_index(p.first, p.second);
  return (size_t) -1;
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::string &marker, double f0, double f2)
{
  return set_ibd( pair_index(i1,i2), marker_index(marker), f0, f2);
}

inline bool
IBD::set_ibd(size_t i, const std::string &marker, double f0, double f2)
{
  return set_ibd( i, marker_index(marker), f0, f2);
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
             size_t m, double f0, double f2)
{
  return set_ibd( pair_index(i1,i2), m, f0, f2);
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::vector<double> &f0, const std::vector<double> &f2)
{
  return set_ibd( pair_index(i1,i2), f0, f2);
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
        const std::string &marker, double f0, double f1, double f2)
{
  return set_ibd( pair_index(i1,i2), marker_index(marker), f0, f1, f2);
}

inline bool
IBD::set_ibd(size_t i, const std::string &marker, double f0, double f1, double f2)
{
  return set_ibd( i, marker_index(marker), f0, f1, f2);
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
        size_t m, double f0, double f1, double f2)
{
  return set_ibd( pair_index(i1,i2), m, f0, f1, f2);
}

inline bool
IBD::set_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::vector<double> &f0, const std::vector<double> &f1, const std::vector<double> &f2)
{
  return set_ibd( pair_index(i1,i2), f0, f1, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::string &marker, double &f0, double &f2) const
{
  return get_ibd( pair_index(i1,i2), marker_index(marker), f0, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             size_t m, double &f0, double &f2) const
{
  return get_ibd( pair_index(i1,i2), m, f0, f2);
}

inline bool
IBD::get_ibd(size_t i, const std::string &marker, double &f0, double &f2) const
{
  return get_ibd( i, marker_index(marker), f0, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             std::vector<double> &f0, std::vector<double> &f2) const
{
  return get_ibd( pair_index(i1,i2), f0, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::string &marker, double &f0, double &f1, double &f2) const
{
  return get_ibd( pair_index(i1,i2), marker_index(marker), f0, f1, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             size_t m, double &f0, double &f1, double &f2) const
{
  return get_ibd( pair_index(i1,i2), m, f0, f1, f2);
}

inline bool
IBD::get_ibd(size_t i, const std::string &marker, double &f0, double &f1, double &f2) const
{
  return get_ibd( i, marker_index(marker), f0, f1, f2);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             std::vector<double> &f0, std::vector<double> &f1, std::vector<double> &f2) const
{
  return get_ibd( pair_index(i1,i2), f0, f1, f2);
}

// Will be taken out of line eventually
inline bool
IBD::get_ibd(size_t i, size_t m, size_t n, double &fn) const
{
  double f0, f1, f2;

  fn = std::numeric_limits<double>::quiet_NaN();

  if(!get_ibd( i, m, f0, f1, f2))
    return false;

  if( SAGE::isnan(f0) || SAGE::isnan(f2))
    return true;

  switch(n)
  {
    case 0: fn = f0;              break;
    case 1: fn = 1.0 - f0 - f2;   break;
    case 2: fn = f2;              break;
    case 3: fn = f1;              break;
  }
  return true;
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             const std::string &marker, size_t n, double &fn) const
{
  return get_ibd( pair_index(i1,i2), marker_index(marker), n, fn);
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             size_t m, size_t n, double &fn) const
{
  return get_ibd( pair_index(i1,i2), m, n, fn);
}

inline bool
IBD::get_ibd(size_t i, const std::string &marker, size_t n, double &fn) const
{
  return get_ibd( i, marker_index(marker), n, fn);
}


// Will be taken out of line eventually
inline bool
IBD::get_ibd(size_t i, size_t n, std::vector<double> &fn) const
{
  std::vector<double> f;

  if(n == 0)
    return get_ibd(i, fn, f);
  else if(n == 2)
    return get_ibd(i, f, fn);
  else if(n == 1 && get_ibd(i, f, fn))
  {
    for(size_t j = 0; j < marker_count(); ++j)
    {
      if(SAGE::isnan(fn[i]) || SAGE::isnan(f[i]))
        fn[i] = std::numeric_limits<double>::quiet_NaN();
      else
        fn[i] = 1.0 - fn[i] - f[i];
    }
    return true;
  }
  else if(n == 3)
  {
    std::vector<double> f2;
    return get_ibd(i, f, fn, f2);
  }
  
  // Catch any errors
  fn.resize(0);
  fn.resize(marker_count(), std::numeric_limits<double>::quiet_NaN());
  return false;
}

inline bool
IBD::get_ibd(const mem_pointer i1, const mem_pointer i2,
             size_t n, std::vector<double> &fn) const
{
  return get_ibd( pair_index(i1,i2), n, fn);
}

inline void
IBD::set_ibd_option(const ibd_option_type& opt)
{
  my_ibd_option = opt;
}

inline
const ibd_option_type&
IBD::get_ibd_option() const
{
  return my_ibd_option;
}

inline size_t
IBD::add_marker(const std::string &name, double dis, gmodel_type m_type)
{
  if( built() )
    return (size_t) -1;

  ibd_marker_info m_info(name, dis, m_type);

  my_markers.push_back(m_info);

  return marker_count() - 1;
}

inline size_t
IBD::marker_count() const
{
  return my_markers.size();
}

inline size_t
IBD::marker_index(const std::string &name) const
{
  for( size_t i = 0; i < marker_count(); ++i )
    if( name == my_markers[i].name )
      return i;

  return (size_t) -1;
}

inline std::string
IBD::marker_name(size_t m) const
{
  return my_markers[m].name;
}

inline const ibd_marker_info&
IBD::get_marker_info(size_t m) const
{
  return my_markers[m];
}

