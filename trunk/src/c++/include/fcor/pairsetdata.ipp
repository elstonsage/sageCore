// ---------------------------------------------------------------------------
// Inline Implementation of PairSetData
// ---------------------------------------------------------------------------

inline const FcorParser*                  
PairSetData::get_parser() const
{
  return my_parser;
}

inline const subtype_vector&
PairSetData::get_subtypes() const
{
  return my_subtypes;
}
/*
inline const maintype_vector&
PairSetData::get_maintypes() const
{
  return my_maintypes;
}
*/
inline const pairset_info_vector&
PairSetData::get_subtype_info() const
{
  return my_subtype_info;
}

inline const pairset_info_vector&
PairSetData::get_maintype_info() const
{
  return my_maintype_info;
}

inline const pairset_vector&
PairSetData::get_subtype_pairset() const
{
  return my_subtype_pairset;
}

inline const pairset_vector&
PairSetData::get_maintype_pairset() const
{
  return my_maintype_pairset;
}

inline
const main_to_sub_map&
PairSetData::get_pairset_group() const
{
  return my_pairset_group;
}

inline bool
PairSetData::is_main_built() const
{
  return my_main_built;
}

// end of PairSetData Inline Implementation

