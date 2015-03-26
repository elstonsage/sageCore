//==========================================================================
//  File:       fped.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                          yjs Nov. 03
//              Major Revision and generalization                gcw Oct. 04
//
//  Notes:      filtered_pedigree build functions.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mped/mp_utilities.h"
#include "fped/fped_obj.h"

namespace SAGE {
namespace FPED {

void
FilteredMultipedigree::construct()
{
  // Build our multipedigree data structures
  build();
  
  // For each pedigree, figure out the source pedigree (either from the members
  // or doing a pedigree find in the source multipedigree), then call construct
  // on the pedigree info class.
  for(pedigree_iterator ped = pedigree_begin(); ped != pedigree_end(); ++ped)
  {
    // Determine the pedigree source

    // First try our members for a source
    const RPED::RefPedigree* source_ped = NULL;
    
    for(member_const_iterator i = ped->member_begin(); !source_ped && i != ped->member_end(); ++i)
    {
      if(i->info().my_member)
      {
        source_ped = i->info().my_member->pedigree();
      }
    }

    // If the member search didn't turn it up, try looking up the name
    if(!source_ped)
    {
      //lint --e{613}  pointer cannot be null 
      source_ped = my_source->pedigree_find(ped->name());

      // If this didn't work, we've got a ped that we can't locate.  Bail out.
      if(!source_ped)
        SAGE_internal_error();
    }
    
    // Construct the pedigree info
    ped->info().construct(*ped, source_ped);
  }
  
  // Finally, sort me!
  sort_into_descent_order();
} 
  
void FilteredMultipedigree::sort_into_descent_order()
{
  for(pedigree_iterator p = pedigree_begin(); p != pedigree_end(); ++p)
  {
    sort_ped_into_descent_order(*p);
  
    SubpedigreeIterator si = p->subpedigree_begin();
    for( ; si != p->subpedigree_end(); ++si )
    {
      sort_sped_into_descent_order(*si);
    }
  }
}

void FilteredMultipedigree::sort_ped_into_descent_order(pedigree_type& p)
{
  for(size_t midx = 0; midx != p.member_count(); ++midx)
  {
    const Member& m = p.member_index(midx);

    // If the member is a nonfounder, check that they come after their parents
    if( m.parent1() )
    {
      // Get the parent indices
      uint p1 = m.parent1()->index();
      uint p2 = m.parent2()->index();
  
      // Choose the larger
      size_t maxp = max(p1,p2);
  
      // If the member's index is smaller than the larger parent, we must reorder
      // them
      if( m.index() < maxp )
      {
        // swap the member with its parent
        p.info().swap_members(midx, maxp);
        p.member_index_swap  (midx, maxp);

        // Back up a step to recheck the parent, which is now at midx.
        --midx;
      }
    }
  }
}

void FilteredMultipedigree::sort_sped_into_descent_order(subpedigree_type& sp)
{
  for(size_t midx = 0; midx != sp.member_count(); ++midx)
  {
    const Member& m = sp.member_index(midx);
    // If the member is a nonfounder, check that they come after their parents
    if( m.parent1() )
    {
      // Get the parent indices
      uint p1 = m.parent1()->subindex();
      uint p2 = m.parent2()->subindex();
  
      // Choose the larger
      size_t maxp = max(p1,p2);
  
      // If the member's index is smaller than the larger parent, we must reorder
      // them
      if( m.subindex() < maxp )
      {
        // swap the member with its parent
        sp.member_index_swap  (midx, maxp);

        // Back up a step to recheck the parent, which is now at midx.
        --midx;
      }
    }
  }
}

// ----------------------------------------------------------------------
  
void
FilteredPedigreeInfo::construct
    (Pedigree& fped, const RPED::RefPedigree* ref_ped)
{
  my_ped_source  = ref_ped;

  //lint --e{613}  pointer should never be null 
  my_ref_pedinfo = &(ref_ped->info());
  
  my_ref_indices.resize(fped.member_count());
  
  std::fill(my_ref_indices.begin(), my_ref_indices.end(), (size_t) -1);

  for( size_t mem = 0; mem < fped.member_count(); ++mem )
  {
    FilteredMemberInfo& fminfo = fped.member_index(mem).info();
    
    // If the member's info is empty, we set it's ref_pedinfo, and skip over the
    // index setting.
    
    if(fminfo.get_source_member() == NULL)
    {
      fminfo.my_ref_pedinfo = my_ref_pedinfo;

      continue;
    }
    
    // If the member's info wasn't empty, we verify that its pedigree info matches
    // ours and set our index vector to match it.
    
    // Verify that it points to the same pedigree that we're working on.  
    // If it's not, that's a big error
    
    if(fminfo.my_ref_pedinfo != my_ref_pedinfo)
      SAGE_internal_error();
    
    // Set the index for this member based upon that in the member's info
    
    my_ref_indices[mem] = fminfo.my_ref_index;
  }
}

void
FilteredPedigreeInfo::swap_members(size_t mem1, size_t mem2)
{
  if( mem1 >= member_count() || mem2 >= member_count() )
    return;

  swap(my_ref_indices[mem1], my_ref_indices[mem2]);
}

void
FilteredPedigreeInfo::dump_ref_indices() const
{
  cout << "FilteredPedigreeInfo::my_ref_indices : " << endl;
  for( size_t i = 0; i < my_ref_indices.size(); ++i )
    cout << my_ref_indices[i] << "  ";
  cout << endl;
}

//
//---------------------------------------------------------------------------------
//

} // end of namespace FPED
} // end of namespace SAGE

