#ifndef SPBASE_H
#include "mped/spbase.h"
#endif

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: subpedigree_base
//============================================================================
//
inline uint
subpedigree_base::index() const
{
    return my_index;
}

inline uint
subpedigree_base::family_count() const
{
    return my_family_index.size();
}

inline uint
subpedigree_base::member_count() const
{
    return my_member_index.size();
}

inline const string&
subpedigree_base::name() const
{
    return my_name;
}

//----------
//
inline subpedigree_base::pedigree_const_pointer
subpedigree_base::pedigree() const
{
    return my_pedigree;
}

inline subpedigree_base::multipedigree_const_pointer
subpedigree_base::multipedigree() const
{
    return my_pedigree->multipedigree();
}

//----------
//
inline const family_base&
subpedigree_base::family_index(uint i) const
{
    return *(my_family_index[i]);
}

inline const member_base&
subpedigree_base::member_index(uint i) const
{
    return *(my_member_index[i]);
}

//----------
//
inline subpedigree_base::family_const_iterator
subpedigree_base::family_begin() const
{
    return family_const_iterator(my_family_index.begin());
}

inline subpedigree_base::family_const_iterator
subpedigree_base::family_end() const
{
    return family_const_iterator(my_family_index.end());
}

inline subpedigree_base::member_const_iterator
subpedigree_base::member_begin() const
{
    return member_const_iterator(my_member_index.begin());
}

inline subpedigree_base::member_const_iterator
subpedigree_base::member_end() const
{
    return member_const_iterator(my_member_index.end());
}

//----------
//
inline subpedigree_base::pedigree_pointer
subpedigree_base::pedigree()
{
    return my_pedigree;
}

inline subpedigree_base::multipedigree_pointer
subpedigree_base::multipedigree()
{
    return my_pedigree->multipedigree();
}

//----------
//
inline family_base&
subpedigree_base::family_index(uint i)
{
    return *(my_family_index[i]);
}

inline member_base&
subpedigree_base::member_index(uint i)
{
    return *(my_member_index[i]);
}

//----------
//
inline subpedigree_base::family_iterator
subpedigree_base::family_begin()
{
    return family_iterator(my_family_index.begin());
}

inline subpedigree_base::family_iterator
subpedigree_base::family_end() 
{
    return family_iterator(my_family_index.end());
}

inline subpedigree_base::member_iterator
subpedigree_base::member_begin()
{
    return member_iterator(my_member_index.begin());
}

inline subpedigree_base::member_iterator
subpedigree_base::member_end()
{
    return member_iterator(my_member_index.end());
}

//----------
//
inline void
subpedigree_base::clear_index_arrays()
{
    my_member_index.clear();
    my_family_index.clear();
}

inline void
subpedigree_base::reset(pedigree_id pid)
{
    my_pedigree = pid;
}

inline void
subpedigree_base::reset(const string& name, pedigree_id pid)
{
    my_name     = name;
    my_pedigree = pid;
}

inline void
subpedigree_base::set_index(uint i)
{
    my_index = i;
}


//============================================================================
//  IMPLEMENTATION: subped_less
//============================================================================
//
inline bool
subped_less::operator ()(subped_id a, subped_id b) const
{
    return a->name() < b->name();
}


}
}

