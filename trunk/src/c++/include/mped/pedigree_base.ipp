#ifndef SPBASE_H
#include "mped/spbase.h"
#endif

namespace SAGE {
namespace MPED {
  
inline uint
pedigree_base::index() const
{
    return my_index;
}

inline uint
pedigree_base::error_count() const
{
    return my_builder.my_errors.size();
}

inline uint
pedigree_base::warning_count() const
{
    return my_builder.my_warnings.size();
}

inline uint
pedigree_base::family_count() const
{
    return my_fcount;
}

inline uint
pedigree_base::member_count() const
{
    return my_mcount;
}

inline uint
pedigree_base::subpedigree_count() const
{
    return my_scount;
}

inline uint
pedigree_base::unconnected_count() const
{
    return my_ucount;
}

inline const string&
pedigree_base::name() const
{
    return my_name;
}

//----------
//
inline pedigree_base::error_iterator
pedigree_base::error_begin() const
{
    return my_builder.my_errors.begin();
}

inline pedigree_base::error_iterator
pedigree_base::error_end() const
{
    return my_builder.my_errors.end();
}

inline pedigree_base::error_iterator
pedigree_base::warning_begin() const
{
    return my_builder.my_warnings.begin();
}

inline pedigree_base::error_iterator
pedigree_base::warning_end() const
{
    return my_builder.my_warnings.end();
}

inline pedigree_base::multipedigree_const_pointer
pedigree_base::multipedigree() const
{
    return my_mped;
}

//----------
//
inline const subpedigree_base&
pedigree_base::subpedigree_index(uint i) const
{
    return *my_subped_index[i];
}

inline const member_base&
pedigree_base::member_index(uint i) const
{
    return *my_member_index[i];
}

inline const family_base&
pedigree_base::family_index(uint i) const
{
    return *my_family_index[i];
}

//----------
//
inline pedigree_base::family_const_pointer
pedigree_base::family_find(const string& p1, const string& p2) const
{
    return lookup_family(p1, p2);
}

inline pedigree_base::family_const_pointer
pedigree_base::family_find(const member_base& p1, const member_base& p2) const
{
    return lookup_family(p1.name(), p2.name());
}

inline pedigree_base::family_const_pointer
pedigree_base::family_find
(const member_const_iterator& p1, const member_const_iterator& p2) const
{
    return lookup_family(p1->name(), p2->name());
}

inline pedigree_base::member_const_pointer
pedigree_base::member_find(const string& m) const
{
    return lookup_member(m);
}

inline pedigree_base::subpedigree_const_pointer
pedigree_base::subpedigree_find(const string& s) const
{
    return lookup_subped(s);
}

//----------
//
inline pedigree_base::family_const_iterator
pedigree_base::family_begin() const
{
    return family_const_iterator(my_family_index.begin());
}

inline pedigree_base::family_const_iterator
pedigree_base::family_begin(const subpedigree_base& s) const
{
    return s.family_begin();
}

inline pedigree_base::family_const_iterator
pedigree_base::family_begin(const subpedigree_const_iterator& s) const
{
    return s->family_begin();
}

inline pedigree_base::family_const_iterator
pedigree_base::family_end() const
{
    return family_const_iterator(my_family_index.end());
}

inline pedigree_base::family_const_iterator
pedigree_base::family_end(const subpedigree_base& s) const
{
    return s.family_end();
}

inline pedigree_base::family_const_iterator
pedigree_base::family_end(const subpedigree_const_iterator& s) const
{
    return s->family_end();
}

inline pedigree_base::member_const_iterator
pedigree_base::member_begin() const
{
    return member_const_iterator(my_member_index.begin());
}

inline pedigree_base::member_const_iterator
pedigree_base::member_begin(const subpedigree_base& s) const
{
    return s.member_begin();
}

inline pedigree_base::member_const_iterator
pedigree_base::member_begin(const subpedigree_const_iterator& s) const
{
    return s->member_begin();
}

inline pedigree_base::member_const_iterator
pedigree_base::member_end() const
{
    return member_const_iterator(my_member_index.end());
}

inline pedigree_base::member_const_iterator
pedigree_base::member_end(const subpedigree_base& s) const
{
    return s.member_end();
}

inline pedigree_base::member_const_iterator
pedigree_base::member_end(const subpedigree_const_iterator& s) const
{
    return s->member_end();
}

inline pedigree_base::subpedigree_const_iterator
pedigree_base::subpedigree_begin() const
{
    return my_subped_index.begin();
}

inline pedigree_base::subpedigree_const_iterator
pedigree_base::subpedigree_end() const
{
    return my_subped_index.end();
}

inline pedigree_base::member_const_iterator
pedigree_base::unconnected_begin() const
{
    return my_unconnecteds->member_begin();
}

inline pedigree_base::member_const_iterator
pedigree_base::unconnected_end() const
{
    return my_unconnecteds->member_end();
}

inline pedigree_base::multipedigree_pointer
pedigree_base::multipedigree()
{
    return my_mped;
}

//----------
//
inline subpedigree_base&
pedigree_base::subpedigree_index(uint i)
{
    return *my_subped_index[i];
}

inline member_base&
pedigree_base::member_index(uint i)
{
    return *my_member_index[i];
}

inline family_base&
pedigree_base::family_index(uint i)
{
    return *my_family_index[i];
}

//----------
//
inline pedigree_base::family_pointer
pedigree_base::family_find(const string& p1, const string& p2)
{
    return lookup_family(p1, p2);
}

inline pedigree_base::family_pointer
pedigree_base::family_find(const member_base& p1, const member_base& p2)
{
    return lookup_family(p1.name(), p2.name());
}

inline pedigree_base::family_pointer
pedigree_base::family_find
(const member_const_iterator& p1, const member_const_iterator& p2)
{
    return lookup_family(p1->name(), p2->name());
}

inline pedigree_base::member_pointer
pedigree_base::member_find(const string& m)
{
    return lookup_member(m);
}

inline pedigree_base::subpedigree_pointer
pedigree_base::subpedigree_find(const string& s)
{
    return lookup_subped(s);
}

//----------
//
inline pedigree_base::family_iterator
pedigree_base::family_begin()
{
    return family_iterator(my_family_index.begin());
}

inline pedigree_base::family_iterator
pedigree_base::family_begin(const subpedigree_base& s)
{
    return (const_cast<subpedigree_base&>(s)).family_begin();
}

inline pedigree_base::family_iterator
pedigree_base::family_begin(const subpedigree_const_iterator& s)
{
    return (const_cast<subpedigree_base&>(*s)).family_begin();
}

inline pedigree_base::family_iterator
pedigree_base::family_end()
{
    return family_iterator(my_family_index.end());
}

inline pedigree_base::family_iterator
pedigree_base::family_end(const subpedigree_base& s)
{
    return (const_cast<subpedigree_base&>(s)).family_end();
}

inline pedigree_base::family_iterator
pedigree_base::family_end(const subpedigree_const_iterator& s)
{
    return (const_cast<subpedigree_base&>(*s)).family_end();
}

inline pedigree_base::member_iterator
pedigree_base::member_begin()
{
    return member_iterator(my_member_index.begin());
}

inline pedigree_base::member_iterator
pedigree_base::member_begin(const subpedigree_base& s)
{
    return (const_cast<subpedigree_base&>(s)).member_begin();
}

inline pedigree_base::member_iterator
pedigree_base::member_begin(const subpedigree_const_iterator& s)
{
    return (const_cast<subpedigree_base&>(*s)).member_begin();
}

inline pedigree_base::member_iterator
pedigree_base::member_end()
{
    return member_iterator(my_member_index.end());
}

inline pedigree_base::member_iterator
pedigree_base::member_end(const subpedigree_base& s)
{
    return (const_cast<subpedigree_base&>(s)).member_end();
}

inline pedigree_base::member_iterator
pedigree_base::member_end(const subpedigree_const_iterator& s)
{
    return (const_cast<subpedigree_base&>(*s)).member_end();
}

inline pedigree_base::subpedigree_iterator
pedigree_base::subpedigree_begin()
{
    return my_subped_index.begin();
}

inline pedigree_base::subpedigree_iterator
pedigree_base::subpedigree_end()
{
    return my_subped_index.end();
}

inline pedigree_base::member_iterator
pedigree_base::unconnected_begin()
{
    return my_unconnecteds->member_begin();
}

inline pedigree_base::member_iterator
pedigree_base::unconnected_end()
{
    return my_unconnecteds->member_end();
}

//----------
//
inline void
pedigree_base::set_index(uint i)
{
    my_index = i;
}


//============================================================================
//  IMPLEMENTATION: ped_less
//============================================================================
//
inline bool
ped_less::operator ()(pedigree_id a, pedigree_id b) const
{
    return a->name() < b->name();
}

}
}

