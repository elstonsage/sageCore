//==========================================================================
//  File:    mcmc_meiosis_map.cpp
//
//  Author:  Geoff Wedig
//           Yeunjoo Song
//
//  History: 0.1 Initial Implementation                         Sept 25 1997 
//           0.2 Revised and moved into own file                Mar  24 1998 
//           0.3 Split into Marker specific                     Apr  07 1998 
//               and non-marker spec. parts                                     
//
//           1.0 Updated to new libraries                        yjs May. 04
//
//  Notes:   Meiosis Mapping - Ordered meioses and individuals for a
//           pedigree to generate likelihood vectors.
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================                                                                               

#include "mcmc/mcmc_meiosis_map.h"

namespace SAGE
{

namespace MCMC
{

McmcMeiosisMap::McmcMeiosisMap(const FPED::Subpedigree& sp, bool x)
  : my_parent_meioses(),
    my_individual_positions(),
    my_mped(sp.multipedigree()),
    my_subped(&sp),
    my_meiosis_count(0),
    my_founder_count(0),
    my_nonfounder_count(0),
    my_x_linked(x)
{
  count_members();

  set_indices();
}

McmcMeiosisMap::McmcMeiosisMap(const McmcMeiosisMap& mm)
  : my_parent_meioses       (mm.my_parent_meioses),
    my_individual_positions (mm.my_individual_positions),
    my_mped                 (mm.my_mped),
    my_subped               (mm.my_subped),
    my_meiosis_count        (mm.my_meiosis_count),
    my_founder_count        (mm.my_founder_count),
    my_nonfounder_count     (mm.my_nonfounder_count),
    my_x_linked             (mm.my_x_linked)
{ }

McmcMeiosisMap&
McmcMeiosisMap::operator=(const McmcMeiosisMap& mm)
{
  if(&mm != this)
  {
    my_parent_meioses       = mm.my_parent_meioses;
    my_individual_positions = mm.my_individual_positions;
  
    my_mped    = mm.my_mped;
    my_subped  = mm.my_subped;
  
    my_meiosis_count    = mm.my_meiosis_count;

    my_founder_count    = mm.my_founder_count;
    my_nonfounder_count = mm.my_nonfounder_count;
  
    my_x_linked = mm.my_x_linked;
  }

  return *this;
}

/// Counts the number of founders and non-founders in the subpedigree
///
void
McmcMeiosisMap::count_members()
{
  for(FPED::MemberConstIterator ind = my_subped->member_begin();
      ind != my_subped->member_end();
      ++ind )
  {
    if(ind->is_founder())
      ++my_founder_count;
    else
      ++my_nonfounder_count;
  }
}

/// Sets the meiosis indices for each non-founder in the pedigree
///
void
McmcMeiosisMap::set_indices()
{
  my_parent_meioses.resize       (my_subped->member_count());
  my_individual_positions.resize (get_nonfounder_count());

  my_meiosis_count = 0;

#if 0
  cout << endl << "mcmc_meiosis_map::set_indices()... " << endl
       << "start : "
       << "nfbits = " << my_meiosis_count << endl;
#endif

  for(FPED::MemberConstIterator i = my_subped->member_begin();
      i != my_subped->member_end(); ++i )
  {
    if( !i->is_founder() )         // If non-founder
    {
      my_individual_positions[my_meiosis_count /2] = i->subindex();

      my_parent_meioses[i->subindex()].moth = my_meiosis_count++;
      my_parent_meioses[i->subindex()].fath = my_meiosis_count++;
    }
  }

#if 0
  cout << "end : "
       << "nfbits = " << my_meiosis_count << endl;

  dump_map();
#endif
}

void
McmcMeiosisMap::dump_map(ostream& o) const
{
  o << "mcmc_meiosis_map dump :" << endl;

  o << "p_meioses :" << endl;
  for(FPED::MemberConstIterator i = my_subped->member_begin();
      i != my_subped->member_end(); ++i )
  {
    o << i->subindex() << " " << i->name() <<  " : "
         << (signed) my_parent_meioses[i->subindex()].moth << ", " << (signed) my_parent_meioses[i->subindex()].fath
         << endl;
  }

  o << endl;
  o << "individual_positions :" << endl;
  for( size_t i = 0; i < my_individual_positions.size(); ++i )
  {
    o << i << " : " << my_individual_positions[i] << endl;
  }

  o << endl;
  o << "nfbits = " << my_meiosis_count << endl;
  o << "founder_count = " << my_founder_count
    << ", nonfounder_count = " << my_nonfounder_count << endl;
  o << "individual count " << get_individual_count() << endl;
  o << "family count "     << get_family_count() << endl;
}

} // end of namespace MCMC

} // end of namespace SAGE

