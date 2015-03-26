// ============================
// marker_likelihood_calculator
// ============================

inline
marker_likelihood_calculator::~marker_likelihood_calculator()
{ }

inline
bool
marker_likelihood_calculator::valid_likelihood() const
{
  for( size_t m = 0; m < my_data->locus_count(); ++m )
  {
    if( !valid_likelihood(m) ) return false;
  }
  return true;
}

inline
bool
marker_likelihood_calculator::valid_likelihood(size_t m) const
{
  if( !my_data->is_valid_locus(m) ) return true;

  return my_graphs[m].is_pattern_valid();
}

inline
double
marker_likelihood_calculator::log_likelihood() const
{
  if( !valid_likelihood() )
  {
    return numeric_limits<double>::quiet_NaN();
  }

  double d = 0.0;

  for( size_t m = 0; m < my_data->locus_count(); ++m )
  {
    d += log_likelihood(m);
  }

  return d;
}

inline
double marker_likelihood_calculator::log_likelihood(size_t m) const
{
  // Make sure this marker matters
  
  if( !my_data->is_valid_locus(m) )
    return 0.0;

  // Get the current bit pattern, and use the other log_likelihood() function
  // to get the value
  
  const bit_field& bits = my_data->get_indicator(m);

#if 0
  cout << "marker_likelihood_calculator::log_likelihood(m = " << m << ") ";
  print_bit_field(cout, bits);
  cout << endl;
#endif

  return log_likelihood(m, bits);
}

inline
double
marker_likelihood_calculator::log_likelihood(size_t m, const bit_field& bits) const
{
#if 0
  cout << "marker_likelihood_calculator::log_likelihood(m = " << m << ", bit = ";
  print_bit_field(cout, bits);
  cout << endl;
#endif

  // Make sure this marker matters
  
  if( !my_data->is_valid_locus(m) )
    return 0.0;

  // If it doesn't, change the graph to the pattern
  
  bool is_valid = my_graphs[m].set_pattern(bits);

  if( !is_valid )
  {
    // If it's not valid, set the cache to NaN

    my_cache_data[m].cache[bits] = numeric_limits<double>::quiet_NaN();

#if 0
  cout << "** pattern not valid.. return qNAN" << endl;
#endif

    return numeric_limits<double>::quiet_NaN();
  }

  // Check if the cache has it.

  hash_iterator h = my_cache_data[m].cache.find(bits);

  if( h != my_cache_data[m].cache.end() )
  {
    ++my_cache_data[m].hit_count;

#if 0
  cout << ".. has it already .. return (*h)->second.get_log()" << endl;
#endif

    return (*h)->second.get_log();
  }

  // We have a valid likelihood.  Let's calculate!
  
  log_double like = my_graphs[m].calculate_likelihood();

  // Cache it and return it.
  
  my_cache_data[m].cache[bits] = like;
  
#if 0
  cout << "return " << like.get_log() << endl;
#endif

  return like.get_log();
}

inline
marker_likelihood_calculator::CachingData::CachingData()
  : hit_count(0)
{
  cache.resize(12323);
}

