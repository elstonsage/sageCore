#ifndef SPBASE_H
#include "mped/spbase.h"
#endif

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: member_base
//============================================================================
//
inline
member_base::member_base(const string& name, const no_info&)
  : my_mpindex((uint) -1),
    my_index((uint) -1),
    my_subindex((uint) -1),
    my_name(name)
{}

//----------
//
inline uint
member_base::mpindex() const
{
    return my_mpindex;
}

inline uint
member_base::index() const
{
    return my_index;
}

inline uint
member_base::subindex() const
{
    return my_subindex;
}

inline uint
member_base::mate_count() const
{
    return my_mate_count;
}

inline uint
member_base::offspring_count() const
{
    return my_offspring_count;
}

inline uint
member_base::sibling_count() const
{
    return (my_origin) ? (my_origin->offspring_count() - 1) : 0;
}

//----------
//

inline SexCode
member_base::get_effective_sex() const
{
    return MPED::get_effective_sex(my_gender);
}

inline SexCode
member_base::get_detailed_sex() const
{
    return my_gender;
}

inline bool
member_base::is_male() const
{
    return MPED::is_male(my_gender);
}

inline bool
member_base::is_female() const
{
    return MPED::is_female(my_gender);
}

inline bool
member_base::is_sex_unknown() const
{
    return MPED::is_sex_unknown(my_gender);
}

inline const string&
member_base::name() const
{
    return my_name;
}

//----------
inline member_base::subpedigree_const_pointer
member_base::subpedigree() const
{
    return my_subped;
}

inline member_base::pedigree_const_pointer
member_base::pedigree() const
{
    return my_subped->pedigree();
}

inline member_base::multipedigree_const_pointer
member_base::multipedigree() const
{
    return my_subped->pedigree()->multipedigree();
}

inline member_base::family_const_pointer
member_base::family() const
{
    return my_origin;
}

inline member_base::member_const_pointer
member_base::parent1() const
{
    return (my_origin) ? my_origin->parent1() : 0;
}

inline member_base::member_const_pointer
member_base::parent2() const
{
    return (my_origin) ? my_origin->parent2() : 0;
}

inline member_base::member_const_pointer
member_base::get_mother() const
{
    return (my_origin) ? my_origin->get_mother() : 0;
}

inline member_base::member_const_pointer
member_base::get_father() const
{
    return (my_origin) ? my_origin->get_father() : 0;
}

//----------
//
inline member_base::mate_const_iterator
member_base::mate_begin() const
{
    return pedigree()->mate_begin(*this);
}

inline member_base::mate_const_iterator
member_base::mate_end() const
{
    return pedigree()->mate_end(*this);
}

inline member_base::offspring_const_iterator
member_base::offspring_begin(const mate_const_iterator& m) const
{
    return pedigree()->offspring_begin(*this, m->mate());
}

inline member_base::offspring_const_iterator
member_base::offspring_begin(const member_const_iterator& m) const
{
    return pedigree()->offspring_begin(*this, *m);
}

inline member_base::offspring_const_iterator
member_base::offspring_begin(const member_base& m) const
{
    return pedigree()->offspring_begin(*this, m);
}

inline member_base::offspring_const_iterator
member_base::offspring_begin(const string& m) const
{
    return pedigree()->offspring_begin(my_name, m);
}

inline member_base::offspring_const_iterator
member_base::offspring_end() const
{
    return pedigree()->offspring_end();
}

inline member_base::parent_const_iterator
member_base::parent_begin() const
{
    return parent_const_iterator(my_origin);
}

inline member_base::parent_const_iterator
member_base::parent_end() const
{
    return parent_const_iterator();
}

inline member_base::progeny_const_iterator
member_base::progeny_begin() const
{
    return pedigree()->progeny_begin(*this);
}

inline member_base::progeny_const_iterator
member_base::progeny_end() const
{
    return pedigree()->progeny_end(*this);
}

inline member_base::sibling_const_iterator
member_base::sibling_begin() const
{
    member_id   data = (my_origin) ? my_origin->offspring() : 0;

    return sibling_const_iterator(const_cast<member_id>(this), data);
}

inline member_base::sibling_const_iterator
member_base::sibling_end() const
{
    return sibling_const_iterator(const_cast<member_id>(this), 0);
}

//----------
inline member_base::subpedigree_pointer
member_base::subpedigree()
{
    return my_subped;
}

inline member_base::pedigree_pointer
member_base::pedigree()
{
    return my_subped->pedigree();
}

inline member_base::multipedigree_pointer
member_base::multipedigree()
{
    return my_subped->pedigree()->multipedigree();
}

inline member_base::family_pointer
member_base::family()
{
    return my_origin;
}

inline member_base::member_pointer
member_base::parent1()
{
    return (my_origin) ? my_origin->parent1() : 0;
}

inline member_base::member_pointer
member_base::parent2()
{
    return (my_origin) ? my_origin->parent2() : 0;
}

inline member_base::member_pointer
member_base::get_mother()
{
    return (my_origin) ? my_origin->get_mother() : 0;
}

inline member_base::member_pointer
member_base::get_father()
{
    return (my_origin) ? my_origin->get_father() : 0;
}

//----------
//
inline member_base::mate_iterator
member_base::mate_begin()
{
    return pedigree()->mate_begin(*this);
}

inline member_base::mate_iterator
member_base::mate_end()
{
    return pedigree()->mate_end(*this);
}

inline member_base::offspring_iterator
member_base::offspring_begin(const mate_const_iterator& m)
{
    return pedigree()->offspring_begin(*this, m->mate());
}

inline member_base::offspring_iterator
member_base::offspring_begin(const member_const_iterator& m)
{
    return pedigree()->offspring_begin(*this, *m);
}

inline member_base::offspring_iterator
member_base::offspring_begin(const member_base& m)
{
    return pedigree()->offspring_begin(*this, m);
}

inline member_base::offspring_iterator
member_base::offspring_begin(const string& m)
{
    return pedigree()->offspring_begin(my_name, m);
}

inline member_base::offspring_iterator
member_base::offspring_end()
{
    return pedigree()->offspring_end();
}

inline member_base::parent_iterator
member_base::parent_begin()
{
    return parent_iterator(my_origin);
}

inline member_base::parent_iterator
member_base::parent_end()
{
    return parent_iterator();
}

inline member_base::progeny_iterator
member_base::progeny_begin()
{
    return pedigree()->progeny_begin(*this);
}

inline member_base::progeny_iterator
member_base::progeny_end()
{
    return pedigree()->progeny_end(*this);
}

inline member_base::sibling_iterator
member_base::sibling_begin()
{
    member_id   data = (my_origin) ? my_origin->offspring() : 0;

    return sibling_iterator(this, data);
}

inline member_base::sibling_iterator
member_base::sibling_end()
{
    return sibling_iterator(this, 0);
}

//----------
//
inline bool member_base::is_founder() const
{
    return parent1() == NULL;
}

inline bool member_base::is_nonfounder() const
{
    return !is_founder();
}

inline bool member_base::is_connected() const
{
  return is_nonfounder() || is_parent();
}

inline bool member_base::is_unconnected() const
{
  return !is_connected();
}

inline bool member_base::is_parent() const
{
  return mate_count() != 0;
}

inline bool member_base::is_non_parent() const
{
  return !is_parent();
}

//----------
//
inline void
member_base::set_sex(SexCode x)
{
    my_gender = x;
}

//----------
//
inline member_id
member_base::siblings() const
{
    return my_siblings;
}

inline uint
member_base::add_child()
{
    return ++my_offspring_count;
}

inline uint
member_base::add_mate()
{
    return ++my_mate_count;
}

inline void
member_base::set_mpindex(uint i)
{
    my_mpindex = i;
}

inline void
member_base::set_index(uint i)
{
    my_index = i;
}

inline void
member_base::set_origin(family_id F)
{
    my_origin = F;
}

inline void
member_base::set_sibling(member_id P)
{
    my_siblings = P;
}

inline void
member_base::set_subped(subped_id SP)
{
    my_subped = SP;
}

inline void
member_base::set_subindex(uint i)
{
    my_subindex = i;
}

//============================================================================
//  IMPLEMENTATION: member_less
//============================================================================
//
inline bool
member_less::operator ()(member_id a, member_id b) const
{
    return a->name() < b->name();
}

}
}

