#ifndef SPBASE_H
#include "mped/spbase.h"
#endif

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: family_base
//============================================================================
//
inline 
family_base::family_base()
  : my_index((uint) -1),
    my_subindex((uint) -1),
    my_subped(0),
    my_parent1(0),
    my_parent2(0),
    my_mother(0),
    my_father(0),
    my_offspring(0),
    my_offspring_count(0)
{}

inline 
family_base::family_base(member_id P1, member_id P2, const no_info&)
  : my_index((uint) -1),
    my_subindex((uint) -1),
    my_subped(0),
    my_parent1(P1),
    my_parent2(P2),
    my_mother(0),
    my_father(0),
    my_offspring(0),
    my_offspring_count(0)
{
    if (P2->name() < P1->name())
    {
        std::swap(my_parent1, my_parent2);
    }
}

//----------
//
inline uint
family_base::index() const
{
    return my_index;
}

inline uint
family_base::subindex() const
{
    return my_subindex;
}

inline uint
family_base::offspring_count() const
{
    return my_offspring_count;
}

inline const string&
family_base::name1() const
{
    return my_parent1->name();
}

inline const string&
family_base::name2() const
{
    return my_parent2->name();
}

inline family_base::subpedigree_const_pointer
family_base::subpedigree() const
{
    return my_subped;
}

inline family_base::pedigree_const_pointer
family_base::pedigree() const
{
    return my_subped->pedigree();
}

inline family_base::multipedigree_const_pointer
family_base::multipedigree() const
{
    return my_subped->pedigree()->multipedigree();
}

inline family_base::member_const_pointer
family_base::parent1() const
{
    return my_parent1;
}

inline family_base::member_const_pointer
family_base::parent2() const
{
    return my_parent2;
}

inline family_base::member_const_pointer
family_base::get_mother() const
{
    return my_mother;
}

inline family_base::member_const_pointer
family_base::get_father() const
{
    return my_father;
}
//----------
//
inline family_base::offspring_const_iterator
family_base::offspring_begin() const
{
    return pedigree()->offspring_begin(*this);
}

inline family_base::offspring_const_iterator
family_base::offspring_end() const
{
    return pedigree()->offspring_end();
}

inline family_base::parent_const_iterator
family_base::parent_begin() const
{
    return parent_const_iterator(const_cast<family_id>(this));
}

inline family_base::parent_const_iterator
family_base::parent_end() const
{
    return parent_const_iterator();
}

//----------
//
inline family_base::subpedigree_pointer
family_base::subpedigree()
{
    return my_subped;
}

inline family_base::pedigree_pointer
family_base::pedigree()
{
    return my_subped->pedigree();
}

inline family_base::multipedigree_pointer
family_base::multipedigree()
{
    return my_subped->pedigree()->multipedigree();
}

inline family_base::member_pointer
family_base::parent1()
{
    return my_parent1;
}

inline family_base::member_pointer
family_base::parent2()
{
    return my_parent2;
}

inline family_base::member_pointer
family_base::get_mother()
{
    return my_mother;
}

inline family_base::member_pointer
family_base::get_father()
{
    return my_father;
}

//----------
//
inline family_base::offspring_iterator
family_base::offspring_begin()
{
    return pedigree()->offspring_begin(*this);
}

inline family_base::offspring_iterator
family_base::offspring_end()
{
    return pedigree()->offspring_end();
}

inline family_base::parent_iterator
family_base::parent_begin()
{
    return parent_iterator(this);
}

inline family_base::parent_iterator
family_base::parent_end()
{
    return parent_iterator();
}

//----------
//
inline member_id
family_base::offspring() const
{
    return my_offspring;
}

inline void
family_base::set_index(uint i)
{
    my_index = i;
}

inline void
family_base::set_subped(subped_id SP)
{
    my_subped = SP;
}

inline void
family_base::set_subindex(uint i)
{
    my_subindex = i;
}

//============================================================================
//  IMPLEMENTATION: family_less
//============================================================================
//
inline bool
family_less::operator ()(family_id a, family_id b) const
{
    return (a->name1() < b->name1()) ||
           (a->name1() == b->name1()  &&  a->name2() < b->name2());
}


}
}

