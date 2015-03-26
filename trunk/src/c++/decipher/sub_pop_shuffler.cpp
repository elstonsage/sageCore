//============================================================================
// File:      sub_pop_shuffler.cpp
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation.                                Feb. 05
//
// Notes:     This file implements a class for generating and storing
//            shuffled sub_populations information for permutation test.
//
// Copyright (c) 2004 R.C. Elston
// All Rights Reserved
//============================================================================

#include "decipher/sub_pop_shuffler.h"

namespace SAGE
{

namespace DECIPHER
{

sub_pop_shuffler::sub_pop_shuffler()
                : my_valid_data(false)
{}

sub_pop_shuffler::sub_pop_shuffler(const vector<const member_em_phenotype_map*>& sub_pops)
                : my_valid_data(false)
{
  set_sub_pops(sub_pops);
}

void
sub_pop_shuffler::set_sub_pops(const vector<const member_em_phenotype_map*>& sub_pops)
{
  my_whole_pop.resize(0);
  my_new_sub_pops.resize(sub_pops.size());

  for( size_t i = 0; i < sub_pops.size(); ++i )
  {
    pool_members(sub_pops[i]->members(), my_new_sub_pops[i]);
  }

  if( my_whole_pop.size() )
    my_valid_data = true;
}

template<typename RandomAccessIter, typename RandomNumberGenerator>
inline void
do_random_shuffle(RandomAccessIter first, RandomAccessIter last,
                  RandomNumberGenerator& rand)
{
  if( first == last )
    return;

  for( RandomAccessIter i = first + 1; i != last; ++i )
    iter_swap(i, first + rand((i - first) + 1));
}

bool
sub_pop_shuffler::do_shuffle(int seed)
{
  if( !my_valid_data )
    return false;

  if( seed )
    my_randomizer.mt.reseed(100);

  //random_shuffle(my_whole_pop.begin(), my_whole_pop.end(), my_randomizer);
  do_random_shuffle(my_whole_pop.begin(), my_whole_pop.end(), my_randomizer);

  size_t last_sub_pop_end = 0;
  for( size_t i = 0; i < my_new_sub_pops.size(); ++i )
  {
    get_new_members(my_new_sub_pops[i], last_sub_pop_end, my_new_sub_pops[i].size());
    last_sub_pop_end += my_new_sub_pops[i].size();
  }

  return true;
}

void
sub_pop_shuffler::dump(ostream& out) const
{
  out << "Whole population:" << endl;

  for( size_t i = 0; i < my_whole_pop.size(); ++i )
  {
    out << "  pedigree  " << my_whole_pop[i]->pedigree()->name()
        << "  member  " << my_whole_pop[i]->name() << endl;
  }
  out << endl;

  out << "New sub populations:" << endl;
  for( size_t s = 0; s < my_new_sub_pops.size(); ++s )
  {
    out << "  Sub population 1:" << s+1 << endl;
    for( size_t i = 0; i < my_new_sub_pops[s].size(); ++i )
    {
      out << "    pedigree  " << my_new_sub_pops[s][i]->pedigree()->name()
          << "  member  " << my_new_sub_pops[s][i]->name() << endl;
    }
    out << endl;
  }
  out << endl;
}

} // end of namespace DECIPHER
} // end of namespace SAGE
