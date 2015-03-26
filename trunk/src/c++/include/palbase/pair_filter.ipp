//////////////////////////////////////////////////////////////////////////
//                Implementation of pair_filter (Inline)                //
//////////////////////////////////////////////////////////////////////////

inline
pair_filter::pair_filter()
{}

inline
pair_filter::~pair_filter()
{}

inline bool
pair_filter::filter_any() const
{
  return (my_markers.size() || my_traits.size());
}

//
//------------------------------------------------------------------
//

inline
pair_filter::filter_trait::filter_trait()
{
  clear();
}

inline
pair_filter::filter_trait::filter_trait(size_t t, double min, double max)
{
  clear();
  this->min = min;
  this->max = max;
  trait = t;
}

inline
pair_filter::filter_trait::filter_trait(size_t t, size_t a, double min, double max)
{
  clear();
  this->min = min;
  this->max = max;
  trait = t;
  if( a < 3 )
  {
    affection[0] = affection[1] = affection[2] = false;
    affection[a] = true;
  }
}

inline
pair_filter::filter_trait::filter_trait(size_t t, const bool* a, double min, double max)
{
  clear();
  this->min = min;
  this->max = max;
  trait = t;
  if( a )
  {
    affection[0] = a[0];
    affection[1] = a[1];
    affection[2] = a[2];
  }
}

inline void
pair_filter::filter_trait::clear()
{
  trait = (size_t)-1;
  affection[0] = affection[1] = affection[2] = true;
  min = NE_INF;
  max = INF;
}

inline bool
pair_filter::filter_trait::operator==(const pair_filter::filter_trait& m) const
{
  return trait == m.trait;
}

inline bool
pair_filter::filter_trait::operator<(const pair_filter::filter_trait& m) const
{
  return trait < m.trait;
}

//
//------------------------------------------------------------------
//

inline
pair_filter::filter_marker::filter_marker()
{
  clear();
}

inline
pair_filter::filter_marker::filter_marker(size_t m, double tolerance)
{
  clear();
  marker = m;
  this->tolerance = tolerance;
}

inline void
pair_filter::filter_marker::clear()
{
  marker = (size_t)-1;
  tolerance = INF;
}

inline bool
pair_filter::filter_marker::operator==(const pair_filter::filter_marker& m) const
{
  return marker == m.marker;
}

inline bool
pair_filter::filter_marker::operator<(const pair_filter::filter_marker& m) const
{
  return marker < m.marker;
}

