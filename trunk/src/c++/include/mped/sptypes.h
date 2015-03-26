#ifndef _SPTYPES_HPP
#define _SPTYPES_HPP

//============================================================================
//  File:       sptypes.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Notes:      This header defines implementation-specific types used to
//              represent pedigree information.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//

#include <vector>
#include <list>
#include <map>
#include <set>

#include "mped/mpfwd.h"

namespace SAGE {
namespace MPED {

//============================================================================
//  PRIVATE SAGE INTERFACE
//
//      This header file contains declarations and definitions that
//      are used for internal representation only.  They should not be
//      used by SAGE client programs, and should be considered undocumented 
//      and subject to change.
//
//      They are left exposed in this manner for two reasons: first, the
//      unreliability of current compilers with respect to namespaces; and
//      second, the heavy inter-dependence of all the types used by the
//      SAGE pedigree representation.
//============================================================================
//
//- The SAGE implementation currently uses a mutable pointer as a unique 
//  identifier.  These identifiers are only used internally.
//
typedef member_base*            member_id;
typedef family_base*            family_id;
typedef subpedigree_base*       subped_id;
typedef pedigree_base*          pedigree_id;
typedef multipedigree_base*     multiped_id;

/** \brief This struct is used to contain non-typed information about the mates of a given member in a pedigree.
  */
class mate_info_base
{
    typedef family_base     family_type;
    typedef member_base     member_type;

  public:
    mate_info_base(member_id M=0, family_id F=0);

    const family_type&  family() const;
    const member_type&  mate() const;
    const string&       name() const;

    family_type&        family();
    member_type&        mate();

  private:
    member_id   the_mate;
    family_id   the_family;
};


/** \brief This struct is used to contain typed information about the mates of a given member in a pedigree.
  */
template <class GI, class FI, class SI, class PI, class MI>
class mate_info : public mate_info_base
{
  public:
    typedef mate_info_base                  base_type;
    typedef member<GI,FI,SI,PI,MI>          member_type;
    typedef SAGE::MPED::family<GI,FI,SI,PI,MI> family_type;
    // Internal note: The above fully-qualified typename is necessary because the parent class
    // defines a function called family(), which collides with the SAGE::MPED::family type.

  public:
    const family_type&  family() const;
    const member_type&  mate() const;

    family_type&        family();
    member_type&        mate();
};



//- We'll use our own assertion macro for reporting excepitons.
//
#define sage_assert(C,X)    if (!(C)) { throw X; }


//- These types are used to create indexable collections of pedigree
//  data objects.
//
typedef std::vector<member_id>      midx_vector;
typedef std::vector<family_id>      fidx_vector;
typedef std::vector<subped_id>      sidx_vector;
typedef std::vector<pedigree_id>    pidx_vector;


//- These container types are used to hold error information
//  in the pedigree.
//
typedef std::list<error_info>       error_list;


//============================================================================
//  STANDARD SET/MAP DECLARATIONS
//============================================================================
//
//  Classes:    member_less, family_less, subped_less
//
//  Purpose:    Functor to determine an ordering of these types.
//
//  Notes:      The ordering of any of these objects is based upon the
//              object's name.
//
//              In order to maximize inlining, these classes' methods are 
//              defined after the definitions of member_base, family_base,
//              and subped_base.
//----------------------------------------------------------------------------
//

/** \internal
  * \brief Thingy
  */
struct member_less
{
    bool    operator ()(member_id a, member_id b) const;
};


/** \internal
  * \brief Thingy
  */
struct family_less
{
    bool    operator ()(family_id a, family_id b) const;
};


/** \internal
  * \brief Thingy
  */
struct subped_less
{
    bool    operator ()(subped_id a, subped_id b) const;
};


/** \internal
  * \brief Thingy
  */
struct ped_less
{
    bool    operator ()(pedigree_id a, pedigree_id b) const;
};


//- These container types actually hold the graph of the pedigree, with its
//  physical and virtual links between members, families, and subpedigrees.
//
typedef std::map<string,pedigree_id>                        pedigree_map;
typedef std::map<pedigree_id,ped_less>                      pedigree_set;

typedef std::set<member_id,member_less>                     member_set;
typedef member_set::const_iterator                          mset_iterator;

typedef std::set<family_id,family_less>                     family_set;
typedef family_set::const_iterator                          fset_iterator;

typedef std::set<subped_id,subped_less>                     subped_set;
typedef subped_set::const_iterator                          sset_iterator;


//- This container and iterator types are used for information about
//  mates.
//
typedef std::multimap<member_id,mate_info_base,member_less> mate_multimap;
typedef mate_multimap::iterator                             mmap_iterator;
typedef mate_multimap::const_iterator                       const_mmap_iterator;
typedef std::pair<const member_id,mate_info_base>           mate_pair;
typedef std::pair<mmap_iterator,mmap_iterator>              mate_range;
typedef std::pair<const_mmap_iterator,const_mmap_iterator>  const_mate_range;

typedef std::set<string>                                    name_set;

} // End namespace MPED
} // End namespace SAGE

#include "mped/sptypes.ipp"

#endif  //- SPTYPES_H
