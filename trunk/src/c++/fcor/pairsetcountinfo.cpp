//****************************************************************************
//* File:      pairsetcountinfo.cpp                                          *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file defines functions to count pair per pedigree *
//*            for each type of pairset of pairsetdata.                      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/pairsetcountinfo.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     PairSetCountInfo                                             ~
// ~                                                                         ~
// ~ Purpose:   Defines functions to count pair per pedigree for each type.  ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of PairSetCountInfo
// ---------------------------------------------------------------------------

PairSetCountInfo::PairSetCountInfo()
{
  my_pedigree_info_built = false;
}

void
PairSetCountInfo::build_pedigree_info(const pairset_type* pairset, size_t ped_count)
{
  my_pedinfo.resize(0);
  my_pedinfo.resize(ped_count, 0);
  
  pairset_type::const_iterator        p = pairset->begin();
  pairset_type::const_iterator pair_end = pairset->end();
  
  for( ; p != pair_end; ++p )
    my_pedinfo[p->member_pair.first->pedigree()->index()] += 1;

  my_pedigree_info_built = true;
}

void
PairSetCountInfo::build_member_info(const pairset_type* pairset, size_t ped_count)
{
  if( !my_pedigree_info_built )
    build_pedigree_info(pairset, ped_count);
  
  my_total_pair_count       = 0;
  my_distinctive_pair_count = 0;

  my_all_member_count    = 0;
  my_first_member_count  = 0;
  my_second_member_count = 0;

  my_first_meminfo.resize(0);
  my_second_meminfo.resize(0);
  my_meminfo.resize(0);

  my_first_meminfo.resize(my_pedinfo.size());
  my_second_meminfo.resize(my_pedinfo.size());
  my_meminfo.resize(my_pedinfo.size());

  pairset_type::const_iterator        p = pairset->begin();
  pairset_type::const_iterator pair_end = pairset->end();
  
  for( ; p != pair_end; ++p )
  {
    //if( p->pair_weight < std::numeric_limits<double>::epsilon() )
    //  continue;

    ++my_total_pair_count;

    my_first_meminfo [p->member_pair.first ->pedigree()->index()].insert(p->member_pair.first);
    my_second_meminfo[p->member_pair.second->pedigree()->index()].insert(p->member_pair.second);

    my_meminfo[p->member_pair.first->pedigree()->index()].insert(p->member_pair.first);
    my_meminfo[p->member_pair.first->pedigree()->index()].insert(p->member_pair.second);
  }

  for( size_t k = 0; k < my_pedinfo.size(); ++k )
  {
    my_first_member_count  += my_first_meminfo[k].size();
    my_second_member_count += my_second_meminfo[k].size();
    my_all_member_count    += my_meminfo[k].size();
  } 

  my_distinctive_pair_count = std::min(my_first_member_count, my_second_member_count);
}  

void
PairSetCountInfo::build_pedigree_info(const pairset_by_pedigree_type& pairset)
{
  my_pedinfo.resize(0);
  my_pedinfo.resize(pairset.size(), 0);
  
  for( size_t p = 0; p < pairset.size(); ++p )
    my_pedinfo[p] = pairset[p].size();

  my_pedigree_info_built = true;
}

void
PairSetCountInfo::build_member_info(const pairset_by_pedigree_type& pairset, bool is_intraclass)
{
  my_total_pair_count       = 0;
  my_distinctive_pair_count = 0;

  my_all_member_count    = 0;
  my_first_member_count  = 0;
  my_second_member_count = 0;

  my_first_meminfo.resize(0);
  my_second_meminfo.resize(0);
  my_meminfo.resize(0);

  my_first_meminfo.resize(pairset.size());
  my_second_meminfo.resize(pairset.size());
  my_meminfo.resize(pairset.size());

  for( size_t p = 0; p < pairset.size(); ++p )
  {
    for( size_t i = 0; i < pairset[p].size(); ++i )
    {
      //if( pairset[p][i].pair_weight < std::numeric_limits<double>::epsilon() )
      //  continue;

      ++my_total_pair_count;

      my_first_meminfo [p].insert(pairset[p][i].member_pair.first);
      my_second_meminfo[p].insert(pairset[p][i].member_pair.second);

      my_meminfo[p].insert(pairset[p][i].member_pair.first);
      my_meminfo[p].insert(pairset[p][i].member_pair.second);
    }

  }

  for( size_t k = 0; k < pairset.size(); ++k )
  {
    my_first_member_count  += my_first_meminfo[k].size();
    my_second_member_count += my_second_meminfo[k].size();
    my_all_member_count    += my_meminfo[k].size();
  } 

  if( is_intraclass )
  {
    my_total_pair_count /= 2;

    vector< vector< set<const pedigree_member_type*> > > family_meminfo;

    family_meminfo.resize(pairset.size());

    for( size_t p = 0; p < my_first_meminfo.size(); ++p )
    {
      set<const pedigree_member_type*>::const_iterator mem_begin = my_first_meminfo[p].begin();
      set<const pedigree_member_type*>::const_iterator mem_end   = my_first_meminfo[p].end();

      if( mem_begin == mem_end )
        continue;

      size_t family_count = (*mem_begin)->pedigree()->family_count();

      family_meminfo[p].resize(family_count);

      for( ; mem_begin != mem_end; ++mem_begin )
      {

        size_t fam_index = (*mem_begin)->family()->index();
        family_meminfo[p][fam_index].insert(*mem_begin);
      }
    }

    size_t distinctive_count = 0;

    for( size_t p = 0; p < family_meminfo.size(); ++p )
    {
      for( size_t f = 0; f < family_meminfo[p].size(); ++f )
      {
        distinctive_count += ( family_meminfo[p][f].size() / 2 );
      }
    } 

    my_distinctive_pair_count = distinctive_count;
  }
  else
  {
    my_distinctive_pair_count = std::min(my_first_member_count, my_second_member_count);
  }

#if 0
  cout << "total_first_mem_count  = " << my_first_member_count << endl;
  cout << "total_second_mem_count = " << my_second_member_count << endl;
  cout << "total_mem_count        = " << my_all_member_count << endl;
  cout << "total_pair_count       = " << my_total_pair_count << endl;
  cout << "distinctive_pair_count = " << my_distinctive_pair_count << endl;
#endif
}  

void
PairSetCountInfo::view_pedinfo(ostream& out) const
{
  for( size_t k = 0; k < my_pedinfo.size(); ++k )
    out << "pedigree " << k << " pair count = " << my_pedinfo[k] << endl;
  out << endl;
}

void
PairSetCountInfo::view_first_meminfo(ostream& out) const
{
  for( size_t k = 0; k < my_first_meminfo.size(); ++k )
  {
    out << "first member count in pedigree " << k << " = " << my_first_meminfo[k].size() << endl;

    set<const pedigree_member_type *>::const_iterator m_begin = my_first_meminfo[k].begin();
    set<const pedigree_member_type *>::const_iterator m_end   = my_first_meminfo[k].end();

    for( ; m_begin != m_end; ++m_begin )
      out << "  " << (*m_begin)->name();
    out << endl;
  }
  out << endl;
}

void
PairSetCountInfo::view_second_meminfo(ostream& out) const
{
  for( size_t k = 0; k < my_second_meminfo.size(); ++k )
  {
    out << "second member count in pedigree " << k << " = " << my_second_meminfo[k].size() << endl;

    set<const pedigree_member_type *>::const_iterator m_begin = my_second_meminfo[k].begin();
    set<const pedigree_member_type *>::const_iterator m_end   = my_second_meminfo[k].end();

    for( ; m_begin != m_end; ++m_begin )
      out << "  " << (*m_begin)->name();
    out << endl;
  }
  out << endl;
}

// end of PairSetCountInfo Implementation

} // end of namespace FCOR
} // end of namespace SAGE

