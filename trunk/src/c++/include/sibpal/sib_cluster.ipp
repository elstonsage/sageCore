//////////////////////////////////////////////////////////////////
//              Implementation of sib_cluster                   //
//////////////////////////////////////////////////////////////////

inline size_t
sib_cluster::valid_sib_count() const
{
  return my_sibs.size();
}
  
inline size_t
sib_cluster::valid_pair_count() const
{
  if( my_valid_fsib_pairs.size() || my_valid_hsib_pairs.size() )
    return my_valid_fsib_pairs.size() + my_valid_hsib_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_fsib_pair_count() const
{
  if( my_valid_fsib_pairs.size() )
    return my_valid_fsib_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_hsib_pair_count() const
{
  if( my_valid_hsib_pairs.size() )
    return my_valid_hsib_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_fsib_mm_pair_count() const
{
  if( my_valid_fsib_mm_pairs.size() )
    return my_valid_fsib_mm_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_fsib_mf_pair_count() const
{
  if( my_valid_fsib_mf_pairs.size() )
    return my_valid_fsib_mf_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_fsib_ff_pair_count() const
{
  if( my_valid_fsib_ff_pairs.size() )
    return my_valid_fsib_ff_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_hsib_mm_pair_count() const
{
  if( my_valid_hsib_mm_pairs.size() )
    return my_valid_hsib_mm_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_hsib_mf_pair_count() const
{
  if( my_valid_hsib_mf_pairs.size() )
    return my_valid_hsib_mf_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline size_t
sib_cluster::valid_hsib_ff_pair_count() const
{
  if( my_valid_hsib_ff_pairs.size() )
    return my_valid_hsib_ff_pairs.size();

  if( my_filter && my_filter->filter_any() )
    return 0;

  return 0;
}

inline sib_pair
sib_cluster::operator[](size_t i) const
{
  if( i < valid_fsib_pair_count() )
    return sib_pair(my_valid_fsib_pairs[i], my_data);

  else if( i < valid_fsib_pair_count() + valid_hsib_pair_count() )
    return sib_pair(my_valid_hsib_pairs[i - valid_fsib_pair_count()], my_data);

  return sib_pair();
}

inline bool
sib_cluster::operator==(const sib_cluster &rhs) const
{
  return (    my_valid_fsib_pairs == rhs.my_valid_fsib_pairs
           && my_valid_hsib_pairs == rhs.my_valid_hsib_pairs );
}

inline bool
sib_cluster::operator!=(const sib_cluster &rhs) const
{
  return !(*this == rhs);
}

inline const sib_set&
sib_cluster::get_valid_sibs() const
{
  return my_sibs;
}

