// ---------------------------------------------------------------------------
// Inline Implementation of ACCal
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// Inline Implementation of VarCovCal
// ---------------------------------------------------------------------------

inline size_t
VarCovCal::find_reltype(const pairset_info_vector& pinfos, string reltype) const
{
  for( size_t i = 0; i < pinfos.size(); ++i )
  {
    if( toUpper(pinfos[i].gname) == reltype || toUpper(pinfos[i].rcode) == reltype )
      return i;
  }

  return (size_t)-1;
}

// end of Inline Implementation
