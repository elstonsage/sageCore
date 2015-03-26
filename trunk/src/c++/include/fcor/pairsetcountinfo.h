#ifndef PAIRSETCOUNTINFO_H
#define PAIRSETCOUNTINFO_H

//****************************************************************************
//* File:      pairsetcountinfo.h                                            *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines functions to count pair per pedigree *
//*            for each type of pairset of pairsetdata.                      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/definitions.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     PairSetCountInfo                                             ~
// ~                                                                         ~
// ~ Purpose:   Defines functions to count pair per pedigree for each type.  ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PairSetCountInfo
{
  public:

    PairSetCountInfo();
   
    void    build_pedigree_info(const pairset_type* set_pointer, size_t ped_count);
    void    build_member_info  (const pairset_type* set_pointer, size_t ped_count);

    void    build_pedigree_info(const pairset_by_pedigree_type& pairset);
    void    build_member_info  (const pairset_by_pedigree_type& pairset, bool is_intra);

    size_t  pair_count(size_t ped)            const;

    size_t  member_count(size_t ped)          const;
 
    size_t  first_member_count (size_t ped)   const;
    size_t  second_member_count(size_t ped)   const;

    size_t  get_all_member_count()            const;
    size_t  get_first_member_count ()         const;
    size_t  get_second_member_count()         const;

    size_t  get_total_pair_count()            const;
    size_t  get_distinctive_pair_count()      const;

    void    view_pedinfo(ostream& out)        const;
    void    view_first_meminfo(ostream& out)  const;
    void    view_second_meminfo(ostream& out) const;    

  private:

    bool                                        my_pedigree_info_built;
    
    size_t                                      my_total_pair_count;
    size_t                                      my_distinctive_pair_count;

    size_t                                      my_all_member_count;
    size_t                                      my_first_member_count;
    size_t                                      my_second_member_count;

    vector<size_t>                              my_pedinfo;

    vector< set<const pedigree_member_type*> >  my_first_meminfo;
    vector< set<const pedigree_member_type*> >  my_second_meminfo;

    vector< set<const pedigree_member_type*> >  my_meminfo;
};
    
#include "fcor/pairsetcountinfo.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
