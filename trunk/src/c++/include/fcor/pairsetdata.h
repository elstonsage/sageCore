#ifndef PAIRSETDATA_H
#define PAIRSETDATA_H

//****************************************************************************
//* File:      pairsetdata.h                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file stores groups of fcor pair types.            *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/parser.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     PairSetData                                                  ~
// ~                                                                         ~
// ~ Purpose:   Stores groups of fcor pairset.                               ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PairSetData
{
  public:

    PairSetData();
    
    void                           build_pairset(const FcorParser&     fp);
    void                           build_pairset(const RelationMatrix& rm,
                                                 const FcorParser&     fp);

    void                           view_pairset(const pairset_vector&      p_vector,
                                                const pairset_info_vector& p_info,
                                                ostream& out)  const;

    const FcorParser*              get_parser()                const;

    const subtype_vector&          get_subtypes()              const;

    const pairset_info_vector&     get_subtype_info()          const;
    const pairset_info_vector&     get_maintype_info()         const;

    const pairset_vector&          get_subtype_pairset()       const;
    const pairset_vector&          get_maintype_pairset()      const;

    const main_to_sub_map&         get_pairset_group()         const;

    bool                           is_main_built()             const;

  protected:
  
    void                           build_main_pairset();

    void                           set_pair_weight(pairset_vector& p_vector);

    const FcorParser*              my_parser;

    subtype_vector                 my_subtypes;

    pairset_info_vector            my_subtype_info;
    pairset_info_vector            my_maintype_info;

    main_to_sub_map                my_pairset_group;

    pairset_vector                 my_subtype_pairset;
    pairset_vector                 my_maintype_pairset;

    bool                           my_main_built;
};
                                 
#include "fcor/pairsetdata.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
