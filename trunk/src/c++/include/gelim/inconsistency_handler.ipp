#ifndef INCONSISTENCY_HANDLER_H
#include "gelim/inconsistency_handler.h"
#endif

namespace SAGE{

// Inlines

inline
inconsistency_handler::inconsistency_handler()
                     : sex_linked_error_exist(false)
{ }

inline
inconsistency_handler::~inconsistency_handler() { }

inline bool
inconsistency_handler::add_error
    (const family_type& family, size_type m, error_type e)
{
  if( e == xy_linked )
    sex_linked_error_exist = true;

  family_error_type& fam = incon_families[&family];

  if(fam.empty())
  {
    build_family(family, fam);

    incon_families_list.push_back(&family);
  }

  for(family_error_type::iterator i = fam.begin(); i != fam.end(); ++i)
  {
//    if(i->second.size() <= m) i->second.resize(m+1, none);

//    if(i->second[m] == none) i->second[m] = e;
    if( e != none && i->second.find(m) == i->second.end() )
    {
      i->second[m] = e;
    }
  }

  return true;
}

inline
bool inconsistency_handler::add_error
   (const family_type& family, ind_id ind, size_type m, error_type e)
{
  if( e == xy_linked )
    sex_linked_error_exist = true;

  family_error_type& fam = incon_families[&family];

  if(fam.empty())
  {
    build_family(family, fam);

    incon_families_list.push_back(&family);
  }

  family_error_type::iterator i = fam.begin();

  // mother
  //
  if( e == xy_linked && i->first == ind)
  {
//    if(i->second.size() <= m) i->second.resize(m+1, none);

//    if(i->second[m] == none) i->second[m] = e;

    if( e != none && i->second.find(m) == i->second.end() )
      i->second[m] = e;
  }
  else if( e != xy_linked )
  {
//    if(i->second.size() <= m) i->second.resize(m+1, none);

//    if(i->second[m] == none) i->second[m] = e;

    if( e != none && i->second.find(m) == i->second.end() )
      i->second[m] = e;
  }

  // father
  //
  ++i;
  
  if( e == xy_linked && i->first == ind)
  {
//  if(i->second.size() <= m) i->second.resize(m+1, none);

//  if(i->second[m] == none) i->second[m] = e;

    if( e != none && i->second.find(m) == i->second.end() )
      i->second[m] = e;
  }
  else if( e != xy_linked )
  {
//  if(i->second.size() <= m) i->second.resize(m+1, none);

//  if(i->second[m] == none) i->second[m] = e;

    if( e != none && i->second.find(m) == i->second.end() )
      i->second[m] = e;
  }  
  
  for(++i ; i != fam.end(); ++i)
  {
    if(i->first == ind)
    {
//      if(i->second.size() <= m) i->second.resize(m+1, none);

//      if(i->second[m] == none) i->second[m] = e;

      if( e != none && i->second.find(m) == i->second.end() )
        i->second[m] = e;
      break;
    }
  }

  return true;
}

inline void inconsistency_handler::mark_info
    (ped_id p, size_type m, bool inconsistent, bool informative)
{
  ped_incon& ped = incon_pedigrees[p];

  if(m >= ped.checked.size())
  {
    ped.checked.resize(m+1);
    ped.informative.resize(m+1);
    ped.inconsistent.resize(m+1);
  }

  ++ped.checked[m];

  if(informative)
  {
    ++ped.informative[m];
  
    if(inconsistent)
      ++ped.inconsistent[m];
  }
}

inline
void inconsistency_handler::clear()
{
  incon_families  = incon_family_map();
  incon_pedigrees = incon_pedigree_map();
}

inline
bool inconsistency_handler::is_sex_linked_error_exist() const
{
  return sex_linked_error_exist;
}

inline
size_t inconsistency_handler::incon_family_count() const
{
  return incon_families.size();
}

inline
const inconsistency_handler::fam_list&
inconsistency_handler::get_incon_family_list() const
{
  return incon_families_list;
}

inline
const inconsistency_handler::family_error_type&
inconsistency_handler::get_incon_family(const family_type& id) const
{
  return incon_families.find(&id)->second;
}

// iteration

inline inconsistency_handler::incon_family_iterator inconsistency_handler::family_begin() const
{
  return incon_families.begin();
}

inline inconsistency_handler::incon_family_iterator inconsistency_handler::family_end()   const
{
  return incon_families.end();
}

inline inconsistency_handler::incon_pedigree_iterator inconsistency_handler::pedigree_begin() const
{
  return incon_pedigrees.begin();
}

inline inconsistency_handler::incon_pedigree_iterator inconsistency_handler::pedigree_end()   const
{
  return incon_pedigrees.end();
}

inline void inconsistency_handler::build_family
    (const family_type& family, family_error_type& fam)
{
  // Get the parents and sort them

  ind_id mother = family.parent1();
  ind_id father = family.parent2();

  if(mother->is_male() || father->is_female())
    std::swap(mother,father);

  fam.push_back(std::make_pair(mother, error_map()));
  fam.push_back(std::make_pair(father, error_map()));

  family_type::offspring_const_iterator p = family.offspring_begin();

  for( ; p != family.offspring_end(); ++p)
  {
    fam.push_back(std::make_pair(&*p, error_map()));
  }
}

}

