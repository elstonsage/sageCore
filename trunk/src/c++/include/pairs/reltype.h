#ifndef RELATIONTYPE_H
#define RELATIONTYPE_H

//****************************************************************************
//* File:      reltype.h                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song & Kevin Jacobs                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     This source file defines a relationtype class to uniquely     *
//*            represent relationships between pairs of individuals          *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>
#include <algorithm>
#include "pairs/reldefs.h"

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     RelationType                                                 ~
// ~                                                                         ~
// ~ Purpose:   Define relationship between two individuals (indi1, indi2).  ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

enum bridge_type          { SELF_BRIDGE=0, MZTWIN_BRIDGE, DZTWIN_BRIDGE,
                            FULL_SIB_BRIDGE, HALF_SIB_BRIDGE, NO_BRIDGE };

enum parental_phase_type  { UNKNOWN_PPHASE=0, PARENT1_PHASE, PARENT2_PHASE   };

class RelationType
{
  public:
    typedef   unsigned short   distance_type;

    RelationType();

    distance_type         distance1()                  const;
    distance_type         distance2()                  const;
    bridge_type           bridge()                     const;
    parental_phase_type   phase1()                     const;
    parental_phase_type   phase2()                     const;
    bool                  has_offspring()              const;

    void                  set_distance1(distance_type d);
    void                  set_distance2(distance_type d);
    void                  set_bridge(bridge_type b);
    void                  set_phase1(parental_phase_type p);
    void                  set_phase2(parental_phase_type p);
    void                  set_has_offspring(bool o);

    void                  swap_phase();

    std::string           str()                      const;

    bool  operator==(const RelationType&)            const;
    bool  operator!=(const RelationType&)            const;
    bool  operator< (const RelationType&)            const;
    bool  operator<=(const RelationType&)            const;
    bool  operator> (const RelationType&)            const;
    bool  operator>=(const RelationType&)            const;

  private:
    distance_type        my_distance1;
    distance_type        my_distance2;
    bridge_type          my_bridge;
    parental_phase_type  my_phase1;
    parental_phase_type  my_phase2;
    bool                 my_has_offspring;

}; // end of class definition

struct CompleteRelationType
{
    CompleteRelationType();

    CompleteRelationType(const RelationType&, const pedigree_member_list&,
                                              const pedigree_member_list&);

    template <class CMP>
    void normalize(const_pedigree_member_pair &p, CMP cmp);

    void swap_phase();

    RelationType                relationship;
    pedigree_member_list        lineage1;
    pedigree_member_list        lineage2;
};

#include "pairs/reltype.ipp"

} // end of namespace SAGE

#endif
