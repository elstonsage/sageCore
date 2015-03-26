// ============================================================================
// Inline Functions
// ============================================================================

namespace SAGE {
namespace RPED {


// ==============
// Locus Iterator
// ==============

inline
genome_description::locus_iterator::locus_iterator() : m(-1), gd(NULL) { }

inline
genome_description::locus_iterator::locus_iterator(const iterator& b)
    : m(b.m), gd(b.gd) { }

inline
genome_description::locus_iterator::reference
genome_description::locus_iterator::operator*() const
{ return locus_type(m, *gd); }

inline
genome_description::locus_iterator&
genome_description::locus_iterator::operator++()
{ ++m; return *this; }

inline
genome_description::locus_iterator
genome_description::locus_iterator::operator++(int)
{ iterator i = *this; ++m; return i; }

inline
genome_description::locus_iterator&
genome_description::locus_iterator::operator--()
{ --m; return *this; }

inline
genome_description::locus_iterator
genome_description::locus_iterator::operator--(int)
{ iterator i = *this; --m; return i; }

inline
bool genome_description::locus_iterator::operator==(const iterator& x) const
{ return gd == x.gd && m == x.m; }

inline
bool genome_description::locus_iterator::operator!=(const iterator& x) const
{ return m != x.m || gd != x.gd; }

inline
genome_description::locus_iterator::locus_iterator
(long _m, genome_description& _gd) : m(_m), gd(&_gd) { }

// ===========
// region_type
// ===========

inline genome_description::region_type::region_type() : m(-1), gd(NULL) { }

inline const std::string& genome_description::region_type::name() const
{ return gd->regions[m].name; }

inline
const std::string& genome_description::region_type::name(const std::string& s)
{ return gd->regions[m].name = s; }

inline genome_description::region_type&
     genome_description::region_type::operator= (const region_type& r)
{ m = r.m; gd = r.gd; return *this; }

inline bool genome_description::region_type::operator==
  (const region_type& r) const
{ return m == r.m && gd == r.gd; }

inline bool genome_description::region_type::operator!=
  (const region_type& r) const
{ return m != r.m || gd != r.gd; }


inline bool genome_description::region_type::valid() const { return gd; }

inline bool genome_description::region_type::is_x_linked() const { return gd->regions[m].x_linked; }

inline 
size_t genome_description::region_type::index() const { return m; }

inline
genome_description::locus_type
genome_description::region_type::locus(size_t l) const
{ return locus_type(gd->regions[m].locus_loc + l, *gd); }

inline
size_t genome_description::region_type::locus_count() const
{
  if(m == -1) return 0;
  
  if(m != (int) gd->regions.size()-1)
    return gd->regions[m + 1].locus_loc - gd->regions[m].locus_loc;
  else
    return gd->locus_count() - gd->regions[m].locus_loc;
}

inline size_t genome_description::region_type::point_count() const
{ if(m != (int) gd->regions.size()-1)
    return gd->regions[m + 1].point_loc - gd->regions[m].point_loc;
  else
    return gd->point_count() - gd->regions[m].point_loc;
}

inline double genome_description::region_type::locus_location(size_t l) const
{ return locus(l).location(); }

inline
double genome_description::region_type::point_location(size_t p) const
{ return gd->points[gd->regions[m].point_loc + p].location; }

inline genome_description::MapTypeShPtr
   genome_description::region_type::map() const
{ return gd->map(); }

inline genome_description::region_type::iterator
genome_description::region_type::begin() const
{ return iterator(gd->regions[m].locus_loc, *gd); }

inline genome_description::region_type::iterator
genome_description::region_type::end() const
{ if(m != (int) gd->regions.size() - 1)
    return iterator(gd->regions[m+1].locus_loc, *gd);
  else
    return iterator(gd->locus_count(), *gd);
}

inline
genome_description::region_type::region_type(long _m, genome_description& v)
    : m(_m), gd(&v) { }

// ==========
// locus_type
// ==========

inline
bool genome_description::locus_type::operator==(const locus_type& l)
{ return m == l.m && &mpp == &l.mpp; }

inline
const std::string& genome_description::locus_type::name() const
{ return locus()->name(); }

inline
const genome_description::locus_info* genome_description::locus_type::locus() const
{ return mpp.loci[m].linfo; }

inline
double genome_description::locus_type::location() const
{ return mpp.loci[m].location; }

inline size_t genome_description::locus_type::region_index() const
{ return mpp.loci[m].region; }

inline genome_description::region_type
    genome_description::locus_type::region() const
{ return region_type(region_index(), mpp); }

inline
double genome_description::locus_type::point_location(long p) const
{
  size_t pt = mpp.loci[m].point_loc + p;
  if((int) region_index() != mpp.points[pt].region)
    return dlimit::infinity();

  return mpp.points[pt].location;
}

inline
double genome_description::locus_type::point_distance(long p) const
{
  size_t pt = mpp.loci[m].point_loc + p;
  
  if((int) region_index() != mpp.points[pt].region)
    return dlimit::infinity();

  return (p>0) ? mpp.points[pt].location - location()
               : location() - mpp.points[pt].location;
}

inline
double genome_description::locus_type::point_theta(long p) const
{
  return mpp.map()->rec_frac(point_distance(p));
}

inline
double genome_description::locus_type::locus_location(long _m) const
{
  if((int) region_index() != mpp.loci[m + _m].region) return dlimit::infinity();

  return mpp.loci[m + _m].location;
}

inline
double genome_description::locus_type::locus_distance(long _m) const
{
  if((int) region_index() != mpp.loci[m + _m].region) return dlimit::infinity();

  return (_m > 0) ? locus_location(_m) - location()
                  : location() - locus_location(_m);
}

inline
double genome_description::locus_type::locus_theta(long _m) const
{
  return mpp.map()->rec_frac(locus_distance(_m));
}

inline
size_t genome_description::locus_type::interval_point_count(long _m) const
{
  size_t rm = region().locus(0).m;

  if((size_t) (m + _m) >= rm + region().locus_count()) return (size_t) -1;

  return mpp.loci[m + _m].point_loc - mpp.loci[m].point_loc;
}

inline size_t genome_description::locus_type::point_prev_count() const
{ return mpp.loci[m].point_loc 
       - mpp.regions[region_index()].point_loc; }

inline size_t genome_description::locus_type::point_next_count() const
{
  if(region_index() != mpp.regions.size() - 1)
    return mpp.regions[region_index() + 1].point_loc
         - mpp.loci[m].point_loc - 1;
  else
    return mpp.points.size() - mpp.loci[m].point_loc - 1;
}

inline size_t genome_description::locus_type::locus_prev_count() const
{ return m - mpp.regions[region_index()].locus_loc; }

inline size_t genome_description::locus_type::locus_next_count() const
{ return region().locus_count() - m - 1; }

inline genome_description::locus_type
    genome_description::locus_type::next_locus() const
{ return locus_type(m+1, mpp); }

inline genome_description::locus_type
    genome_description::locus_type::prev_locus() const
{ return locus_type(m-1, mpp); }

inline size_t genome_description::locus_type::point_index() const
{ return point_prev_count(); }

inline size_t genome_description::locus_type::locus_index() const
{ return locus_prev_count(); }

inline
genome_description::locus_type::locus_type(long _m, genome_description& _v)
    : m(_m), mpp(_v) { }

// ==================
// genome_description
// ==================
    
inline
bool genome_description::freeze()
{
  if(!built()) return false;

  return my_frozen = true;
}

inline bool genome_description::built()  const { return my_built; }
inline bool genome_description::frozen() const { return my_frozen; }

inline
bool genome_description::set_mapping_function(MapTypeShPtr mt)
{ 
  if(mt == my_map_func)
    return false;

  my_map_func = mt; 
  return true;  
}

inline
genome_description::MapTypeShPtr genome_description::map() const
{ return my_map_func; }

inline
bool genome_description::set_interval_count(long c)
{ if(c < 1) c = 2; interval = c; method = true; return true; }

inline
bool genome_description::set_scan_distance(double d, bool _absolute)
{ if(d <= 0) d = 2.0; distance = d; absolute = _absolute; method = false;
  return true;
}
  
inline
long genome_description::point_count() const
{ return points.size(); }

inline
long genome_description::locus_count() const
{ return loci.size(); }

inline
long genome_description::region_count() const
{ return regions.size(); }

inline genome_description::locus_type genome_description::locus(size_t m)
{ return locus_type(m, *this); }

inline genome_description::region_type genome_description::region(size_t m)
{ return region_type(m, *this); }

inline genome_description::region_type genome_description::region
  (const std::string& s)
{
  for(int i = 0; i < region_count(); ++i)
  {
    if(toUpper(region_name(i)) == toUpper(s)) return region(i);
  }
  return region_type();
}

inline const std::string& genome_description::region_name(size_t m) const
{ return regions[m].name; }

inline const std::string&
   genome_description::region_name(size_t m, const std::string& s)
{ return regions[m].name = s; }

inline genome_description::iterator genome_description::begin()
{ return iterator(0, *this); }

inline genome_description::iterator genome_description::end()
{ return iterator(locus_count(), *this); }

// ================================
// genome_description::point_struct
// ================================

inline genome_description::point_struct::point_struct()
{  
  location = numeric_limits<double>::quiet_NaN();
  locus = -1;   
}

} // End namespace RPED
} // End namespace SAGE
