// ================================
// Inline functions of ibd_analysis
// ================================

inline bool
ibd_analysis::valid() const
{
  return my_valid;
}

inline double
ibd_analysis::prob_share(mem_pointer p1, mem_pointer p2, long n) const
{
  if( !valid() )
    return numeric_limits<double>::quiet_NaN();
  
  size_t x = (signed) find_id(p1);
  size_t y = (signed) find_id(p2);

  if( x == (size_t) -1 || y == (size_t) -1 )
    return numeric_limits<double>::quiet_NaN();

  // Make p temporary
  {
    size_t p = max(x,y);
  
    x = min(x,y);
    y = p;
  }

  double d;

  switch(n)
  {
    case 0 : return my_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];

    case 1 : d = 1.0 - my_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x]
                     - my_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
             return (d < 0) ? 0.0 : d;

    case 2 : return my_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];

    default: return numeric_limits<double>::quiet_NaN();
  }
}

inline double
ibd_analysis::prob_share(long x, long y, long n) const
{
  if( !valid() )
    return numeric_limits<double>::quiet_NaN();
  
  long p = max(x,y);
  
  x = min(x,y);
  y = p;

  // If x or y is out of bounds, we return a nan.  This can easily be checked
  // by checking if the unsigned x or y is greater than the member count
  if( (size_t) x >= my_meiosis_map.get_subpedigree()->member_count() ||
      (size_t) y >= my_meiosis_map.get_subpedigree()->member_count()    )
    return numeric_limits<double>::quiet_NaN();

  double d;

  switch(n)
  {
    case 0 : return my_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];

    case 1 : d = 1 - my_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x]
                   - my_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
             return (d < 0) ? 0.0 : d;

    case 2 : return my_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];

    default: return numeric_limits<double>::quiet_NaN();
  }
}

inline double
ibd_analysis::sib_prob_share(mem_pointer p1, mem_pointer p2, long n) const
{
  if( !valid() )
    return numeric_limits<double>::quiet_NaN();
  
  size_t x = (signed) find_id(p1);
  size_t y = (signed) find_id(p2);

  if( x == (size_t) -1 || y == (size_t) -1 )
    return numeric_limits<double>::quiet_NaN();

  if( !is_sib(my_meiosis_map, x, y) && !is_hsib(my_meiosis_map, x, y))
    return numeric_limits<double>::quiet_NaN();

  // Make p temporary
  {
    size_t p = max(x,y);
  
    x = min(x,y);
    y = p;
  }

  double d;

  switch(n)
  {
    case 0 : return my_sib_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];

    case 1 : d = 1.0 - my_sib_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x]
                     - my_sib_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
             return (d < 0) ? 0.0 : d;

    case 2 : return my_sib_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];

    default: return numeric_limits<double>::quiet_NaN();
  }
}

inline double
ibd_analysis::sib_prob_share(long x, long y, long n) const
{
  if( !valid() )
    return numeric_limits<double>::quiet_NaN();

  if( !is_sib(my_meiosis_map, x, y) && !is_hsib(my_meiosis_map, x, y) )
    return numeric_limits<double>::quiet_NaN();
  
  long p = max(x,y);
  
  x = min(x,y);
  y = p;

  // If x or y is out of bounds, we return a nan.  This can easily be checked
  // by checking if the unsigned x or y is greater than the member count
  if( (size_t) x >= my_meiosis_map.get_subpedigree()->member_count() ||
      (size_t) y >= my_meiosis_map.get_subpedigree()->member_count()    )
    return numeric_limits<double>::quiet_NaN();

  double d;

  switch(n)
  {
    case 0 : return my_sib_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];

    case 1 : d = 1 - my_sib_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x]
                   - my_sib_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
             return (d < 0) ? 0.0 : d;

    case 2 : return my_sib_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];

    default: return numeric_limits<double>::quiet_NaN();
  }
}
/*
inline double
ibd_analysis::get_mean_share_covariance(mem_pointer p1, mem_pointer p2, mem_pointer q1, mem_pointer q2) const
{
  if( !valid() )
    return numeric_limits<double>::quiet_NaN();
  
  size_t x1 = (signed) find_id(p1);
  size_t y1 = (signed) find_id(p2);
  size_t x2 = (signed) find_id(q1);
  size_t y2 = (signed) find_id(q2);

  if(    x1 == (size_t) -1 || y1 == (size_t) -1
      || x2 == (size_t) -1 || y2 == (size_t) -1 )
    return numeric_limits<double>::quiet_NaN();

  size_t member_count = my_meiosis_map.get_subpedigree()->member_count();

  size_t r_off_set = (min(x1,y1)+1) * (min(x1,y1)+2) / 2;
  size_t c_off_set = (min(x2,y2)+1) * (min(x2,y2)+2) / 2;

  size_t r = min(x1,y1) * member_count + max(x1,y1) - r_off_set;
  size_t c = min(x2,y2) * member_count + max(x2,y2) - c_off_set;

  return my_mean_share_info.covariance(r, c);
}
*/
inline void
ibd_analysis::get_lvec_probability(vector<double>& lvec_prob) const
{
  if( !valid() )
    return;

  lvec_prob.resize(0);

  for( size_t i = 0; i < my_ibd_state.size(); ++i )
    lvec_prob.push_back(my_ibd_state[i].first);

  return;
}

inline void
ibd_analysis::get_pair_ibd_state(mem_pointer p1, mem_pointer p2, vector<size_t>& pair_ibd_state) const
{
  if( !valid() )
    return;

  size_t x = (signed) find_id(p1);
  size_t y = (signed) find_id(p2);

  if( x == (size_t) -1 || y == (size_t) -1 )
    return;

  size_t member_count = my_meiosis_map.get_subpedigree()->member_count();

  size_t r_off_set = (min(x,y)+1) * (min(x,y)+2) / 2;

  size_t r = min(x,y) * member_count + max(x,y) - r_off_set;

  pair_ibd_state.resize(0);

  for( size_t i = 0; i < my_ibd_state.size(); ++i )
    pair_ibd_state.push_back(my_ibd_state[i].second[r]);

  return;
}

inline const a_marker_ibd_state&
ibd_analysis::get_ibd_state() const
{
  return my_ibd_state;
}

inline size_t
ibd_analysis::find_id(mem_pointer i) const
{
  return valid() ? i->subindex() : (size_t) -1;
}

inline bool
ibd_analysis::is_sib(const meiosis_map& mm, long j, long k) const
{
  return is_fsib(mm.member(j), mm.member(k));
}

inline bool
ibd_analysis::is_hsib(const meiosis_map& mm, long j, long k) const
{
  return SAGE::is_hsib(mm.member(j), mm.member(k));
}

inline bool
ibd_analysis::is_maternal_hsib(const meiosis_map& mm, long j, long k) const
{
  return SAGE::is_maternal_hsib(mm.member(j), mm.member(k));
}

inline bool
ibd_analysis::is_paternal_hsib(const meiosis_map& mm, long j, long k) const
{
  return SAGE::is_paternal_hsib(mm.member(j), mm.member(k));
}

inline bool
ibd_analysis::is_brother_brother(const meiosis_map& mm, long j, long k) const
{
  return SAGE::is_brother_brother(mm.member(j), mm.member(k));
}

inline double&
ibd_analysis::my_prob_share(long x, long y, long n)
{
  if(n == 0) return my_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
  else       return my_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];
}

inline double&
ibd_analysis::my_sib_prob_share(long x, long y, long n)
{
  if(n == 0) return my_sib_ibd_sharing[x * my_meiosis_map.get_subpedigree()->member_count() + y];
  else       return my_sib_ibd_sharing[y * my_meiosis_map.get_subpedigree()->member_count() + x];
}
