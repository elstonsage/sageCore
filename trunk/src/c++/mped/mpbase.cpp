//============================================================================
//  File:       mpbase.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <mpfwd.h>
    #pragma hdrstop
#endif

#include "mped/mpbase.h"

namespace SAGE {
namespace MPED {

multipedigree_base::multipedigree_base()
{}


multipedigree_base::~multipedigree_base()
{
    if (pedigree_count() > 0)
    {
        clear();
    }
}


//----------
//
multipedigree_base::member_const_pointer
multipedigree_base::member_find(const string& ped, const string& name) const
{
    pedigree_const_pointer  pp = pedigree_find(ped);

    return (pp) ? pp->member_find(name) : 0;
}


multipedigree_base::pedigree_const_pointer
multipedigree_base::pedigree_find(const string& ped) const
{
    pedigree_map::const_iterator    iter = my_pedmap.find(ped);

    return (iter != my_pedmap.end()) ? iter->second : 0;
}


//----------
//
void
multipedigree_base::pedigree_index_swap(uint i, uint j)
{
    std::swap(my_ped_index[i], my_ped_index[j]);
    my_ped_index[i]->set_index(i);
    my_ped_index[j]->set_index(j);
}


//----------
//
multipedigree_base::member_pointer
multipedigree_base::member_find(const string& ped, const string& name)
{
    pedigree_pointer    pp = pedigree_find(ped);

    return (pp) ? pp->member_find(name) : 0;
}


multipedigree_base::pedigree_pointer
multipedigree_base::pedigree_find(const string& ped)
{
    pedigree_map::iterator  iter = my_pedmap.find(ped);

    return (iter != my_pedmap.end()) ? iter->second : 0;
}


//----------
//
void
multipedigree_base::member_index_swap(uint i, uint j)
{
    std::swap(my_mem_mpindex[i], my_mem_mpindex[j]);
    my_mem_mpindex[i]->set_mpindex(i);
    my_mem_mpindex[j]->set_mpindex(j);
}


//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::add_lineage
(const string& ped, const string& child, const string& parent)
{
    update_cache(ped)->add_lineage(child, parent);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_lineage
(const pedigree_iterator& ped, const string& child, const string& parent)
{
    update_cache(ped)->add_lineage(child, parent);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_lineage
(const string& ped, const string& child, const string& parent1, const string& parent2)
{
    update_cache(ped)->add_lineage(child, parent1, parent2);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_lineage
(const pedigree_iterator& ped, const string& child, const string& parent1, const string& parent2)
{
    update_cache(ped)->add_lineage(child, parent1, parent2);
    return my_last_iter;
}


//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::add_marriage
(const string& ped, const string& spouse1, const string& spouse2)
{
    update_cache(ped)->add_marriage(spouse1, spouse2);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_marriage
(const pedigree_iterator& ped, const string& spouse1, const string& spouse2)
{
    update_cache(ped)->add_marriage(spouse1, spouse2);
    return my_last_iter;
}

//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::add_member
(const string& ped, const string& name, SexCode s)
{
    update_cache(ped)->add_member(name, s);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_member
(const pedigree_iterator& ped, const string& name, SexCode s)
{
    update_cache(ped)->add_member(name, s);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_member
(const string& ped, const string& name, SexCode s, const void* info)
{
    update_cache(ped)->add_member(name, s, info);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_member
(const pedigree_iterator& ped, const string& name, SexCode s, const void* info)
{
    update_cache(ped)->add_member(name, s, info);
    return my_last_iter;
}

//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::add_sibship
(const string& ped, const string& sib1, const string& sib2)
{
    update_cache(ped)->add_sibship(sib1, sib2);
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::add_sibship
(const pedigree_iterator& ped, const string& sib1, const string& sib2)
{
    update_cache(ped)->add_sibship(sib1, sib2);
    return my_last_iter;
}

//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::build(const string& ped)
{
    update_cache(ped)->build();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::build(const pedigree_iterator& ped)
{
    update_cache(ped)->build();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::clear(const string& ped)
{
    update_cache(ped)->clear();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::clear(const pedigree_iterator& ped)
{
    update_cache(ped)->clear();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::flush(const string& ped)
{
    update_cache(ped)->flush_build_info();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::flush(const pedigree_iterator& ped)
{
    update_cache(ped)->flush_build_info();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::freeze(const string& ped)
{
    update_cache(ped)->freeze();
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::freeze(const pedigree_iterator& ped)
{
    update_cache(ped)->freeze();
    return my_last_iter;
}


//----------------------------------------------------------------------------
//
pedigree_id
multipedigree_base::allocate_pedigree(const string& name)
{
    return new pedigree_base(name);
}


void
multipedigree_base::deallocate_pedigree(pedigree_id P)
{
    delete P;
}


void
multipedigree_base::build()
{
    pedigree_map::iterator  pf = my_pedmap.begin();
    pedigree_map::iterator  pl = my_pedmap.end();

    for (;  pf != pl;  ++pf)
    {
        pf->second->build();
    }

    // Assign indices for each member, now that the pedigrees are done.

    assign_member_indices();
}


void
multipedigree_base::clear()
{
    pedigree_map::iterator  pf = my_pedmap.begin();
    pedigree_map::iterator  pl = my_pedmap.end();

    for (;  pf != pl;  ++pf)
    {
        deallocate_pedigree(pf->second);
        pf->second = 0;
    }
    my_pedmap.clear();
    my_ped_index.clear();
}


void
multipedigree_base::flush()
{
    pedigree_map::iterator  pf = my_pedmap.begin();
    pedigree_map::iterator  pl = my_pedmap.end();

    for (;  pf != pl;  ++pf)
    {
        pf->second->flush_build_info();
    }
}


void
multipedigree_base::freeze()
{
    pedigree_map::iterator  pf = my_pedmap.begin();
    pedigree_map::iterator  pl = my_pedmap.end();

    for (;  pf != pl;  ++pf)
    {
        pf->second->freeze();
    }
}


//----------
//
multipedigree_base::pedigree_iterator
multipedigree_base::update_cache(const string& ped_name)
{
    if (ped_name != my_last_name)
    {
        pedigree_map::iterator  iter = my_pedmap.find(ped_name);

        if (iter == my_pedmap.end())
        {
            pedigree_id     id = allocate_pedigree(ped_name);
            insert_pair     vp = my_pedmap.insert( value_pair(ped_name,id) );

            iter = vp.first;
            update_indices(id);
        }
        
        my_last_iter = my_ped_index.begin() + iter->second->index();
        my_last_name = ped_name;
    }
    return my_last_iter;
}


multipedigree_base::pedigree_iterator
multipedigree_base::update_cache(const pedigree_iterator& ped)
{
    if (ped != my_last_iter)
    {
        my_last_iter = ped;
        my_last_name = ped->name();
    }
    return my_last_iter;
}


void
multipedigree_base::update_indices(pedigree_id P)
{
    my_ped_index.push_back(P);
    my_ped_index.back()->set_index( my_ped_index.size() - 1 );
}

void
multipedigree_base::assign_member_indices()
{
  // Get the current member count.  New members are initially assigned ids
  // after the current batch.
  
  uint prior_member_count = member_count();
  
  // Determine the number of members by iterating through the pedigrees and
  // summing
  
  uint total_members = 0;

  for(pedigree_const_iterator i = pedigree_begin(); i != pedigree_end(); ++i)
  {
    total_members += i->member_count();
  }

  // Resize the member index vector
  
  my_mem_mpindex.resize(total_members);
  
  // Go through the pedigrees, looking for members that aren't assigned ids and
  // assign them one

  uint current_id = prior_member_count;
  
  for(pedigree_iterator i = pedigree_begin(); i != pedigree_end(); ++i)
  {
    for(member_iterator j = i->member_begin(); j != i->member_end(); ++j)
    {
      if(j->mpindex() == (uint) -1)
      {
        my_mem_mpindex[current_id] = &*j;
        j->set_mpindex(current_id);
        ++current_id;
      }
    }
  }
  
  // Verify that all ids have been assigned.  If not, there's something seriously
  // wrong.
  
  assert(current_id == total_members);
}

} // End namespace MPED
} // End namespace SAGE
