// ---------------------------------------------------------------------------
// Inline Implementation of PairSetCountInfo
// ---------------------------------------------------------------------------

inline
size_t
PairSetCountInfo::pair_count(size_t ped) const
{
  return my_pedinfo[ped];
}

inline
size_t
PairSetCountInfo::first_member_count(size_t ped) const
{
  return my_first_meminfo[ped].size();
}

inline
size_t
PairSetCountInfo::second_member_count(size_t ped) const
{
  return my_second_meminfo[ped].size();
}

inline
size_t
PairSetCountInfo::member_count(size_t ped) const
{
  return my_meminfo[ped].size();
}

inline
size_t
PairSetCountInfo::get_all_member_count() const
{
  return my_all_member_count;
}

inline
size_t
PairSetCountInfo::get_first_member_count() const
{
  return my_first_member_count;
}

inline
size_t
PairSetCountInfo::get_second_member_count() const
{
  return my_second_member_count;
}

inline
size_t
PairSetCountInfo::get_total_pair_count() const
{
  return my_total_pair_count;
}

inline
size_t
PairSetCountInfo::get_distinctive_pair_count() const
{
  return my_distinctive_pair_count;
}

// end of PairSetCountInfo Implementation

