//============================================================================
// IMPLEMENTATION:  sub_pop_shuffler
//============================================================================
//

inline const vector<member>&
sub_pop_shuffler::get_whole_pop() const
{
  return my_whole_pop;
}

inline const sub_pop_shuffler::sub_pops&
sub_pop_shuffler::get_new_sub_pops() const
{
  return my_new_sub_pops;
}

inline void
sub_pop_shuffler::pool_members(const set<member, member_order<member> >& mems,
                               vector<member>& new_mems)
{
  set<member, member_order<member> >::const_iterator mi = mems.begin();
  for( ; mi != mems.end(); ++mi )
  {
    my_whole_pop.push_back(*mi);
    new_mems.push_back(*mi);
  }

  return;
}

inline void
sub_pop_shuffler::get_new_members(vector<member>& new_sub, size_t start, size_t count)
{
  for(size_t i = 0; i < count; ++i )
    new_sub[i] = my_whole_pop[i + start];

  return;
}

