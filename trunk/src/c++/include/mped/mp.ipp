//============================================================================
//  File:       mp.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//
//============================================================================
//  IMPLEMENTATION: MULTIPEDIGREE
//============================================================================
//

namespace SAGE {
namespace MPED {

template <class GI, class FI, class SI, class PI, class MI>
multipedigree<GI,FI,SI,PI,MI>::multipedigree()
{}

template <class GI, class FI, class SI, class PI, class MI>
multipedigree<GI,FI,SI,PI,MI>::~multipedigree()
{
    clear();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline uint 
multipedigree<GI,FI,SI,PI,MI>::pedigree_count() const
{
    return multipedigree_base::pedigree_count();
}

template <class GI, class FI, class SI, class PI, class MI> inline const MI& 
multipedigree<GI,FI,SI,PI,MI>::info() const
{
    return my_info;
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline uint 
multipedigree<GI,FI,SI,PI,MI>::member_count() const
{
    return multipedigree_base::member_count();
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
const typename multipedigree<GI,FI,SI,PI,MI>::pedigree_type&
multipedigree<GI,FI,SI,PI,MI>::pedigree_index(uint i) const
{
    return reinterpret_cast<const pedigree_type&>(base::pedigree_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_const_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_begin() const
{
    return pedigree_const_iterator(base::pedigree_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_const_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_end() const
{
    return pedigree_const_iterator(base::pedigree_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_const_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_last() const
{
    return pedigree_const_iterator(base::pedigree_last());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
const typename multipedigree<GI,FI,SI,PI,MI>::member_type&
multipedigree<GI,FI,SI,PI,MI>::member_index(uint i) const
{
    return reinterpret_cast<const member_type&>(base::member_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::member_const_pointer
multipedigree<GI,FI,SI,PI,MI>::member_find
(const string& ped, const string& name) const
{
    return reinterpret_cast<member_const_pointer>(base::member_find(ped,name));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_const_pointer
multipedigree<GI,FI,SI,PI,MI>::pedigree_find(const string& name) const
{
    return reinterpret_cast<pedigree_const_pointer>(base::pedigree_find(name));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_type&
multipedigree<GI,FI,SI,PI,MI>::pedigree_index(uint i)
{
    return reinterpret_cast<pedigree_type&>(base::pedigree_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::pedigree_index_swap(uint i, uint j)
{
    base::pedigree_index_swap(i, j);
}

template <class GI, class FI, class SI, class PI, class MI> inline MI& 
multipedigree<GI,FI,SI,PI,MI>::info()
{
    return my_info;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_begin()
{
    return pedigree_iterator(base::pedigree_begin());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_end()
{
    return pedigree_iterator(base::pedigree_end());
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::pedigree_last()
{
    return pedigree_iterator(base::pedigree_last());
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::member_type&
multipedigree<GI,FI,SI,PI,MI>::member_index(uint i)
{
    return reinterpret_cast<member_type&>(base::member_index(i));
}

template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::member_index_swap(uint i, uint j)
{
    base::member_index_swap(i, j);
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::member_pointer
multipedigree<GI,FI,SI,PI,MI>::member_find
(const string& ped, const string& name)
{
    return reinterpret_cast<member_pointer>(base::member_find(ped,name));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_pointer
multipedigree<GI,FI,SI,PI,MI>::pedigree_find(const string& name)
{
    return reinterpret_cast<pedigree_pointer>(base::pedigree_find(name));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_lineage
(const string& ped, const string& child, const string& parent)
{
    return pedigree_iterator(base::add_lineage(ped, child, parent));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_lineage
(const pedigree_iterator& ped, const string& child, const string& parent)
{
    return pedigree_iterator(base::add_lineage(ped, child, parent));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_lineage
(const string& ped, const string& child, const string& parent1, const string& parent2)
{
    return pedigree_iterator(base::add_lineage(ped, child, parent1, parent2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_lineage
(const pedigree_iterator& ped, const string& child, 
 const string& parent1, const string& parent2)
{
    return pedigree_iterator(base::add_lineage(ped, child, parent1, parent2));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_marriage
(const string& ped, const string& spouse1, const string& spouse2)
{
    return pedigree_iterator(base::add_marriage(ped, spouse1, spouse2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_marriage
(const pedigree_iterator& ped, const string& spouse1, const string& spouse2)
{
    return pedigree_iterator(base::add_marriage(ped, spouse1, spouse2));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_member
(const string& ped, const string& name, SexCode s)
{
    return pedigree_iterator(base::add_member(ped, name, s));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_member
(const pedigree_iterator& ped, const string& name, SexCode s)
{
    return pedigree_iterator(base::add_member(ped, name, s));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_member
(const string& ped, const string& name, SexCode s, const geninfo_type& info)
{
    return pedigree_iterator(base::add_member(ped, name, s, &info));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_member
(const pedigree_iterator& ped, const string& name, SexCode s, const geninfo_type& info)
{
    return pedigree_iterator(base::add_member(ped, name, s, &info));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_sibship
(const string& ped, const string& sib1, const string& sib2)
{
    return pedigree_iterator(base::add_sibship(ped, sib1, sib2));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::add_sibship
(const pedigree_iterator& ped, const string& sib1, const string& sib2)
{
    return pedigree_iterator(base::add_sibship(ped, sib1, sib2));
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::build(const string& ped)
{
    return pedigree_iterator(base::build(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::build(const pedigree_iterator& ped)
{
    return pedigree_iterator(base::build(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::clear(const string& ped)
{
    return pedigree_iterator(base::clear(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::clear(const pedigree_iterator& ped)
{
    return pedigree_iterator(base::clear(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::flush(const string& ped)
{
    return pedigree_iterator(base::flush(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::flush(const pedigree_iterator& ped)
{
    return pedigree_iterator(base::flush(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::freeze(const string& ped)
{
    return pedigree_iterator(base::freeze(ped));
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename multipedigree<GI,FI,SI,PI,MI>::pedigree_iterator
multipedigree<GI,FI,SI,PI,MI>::freeze(const pedigree_iterator& ped)
{
    return pedigree_iterator(base::freeze(ped));
}


//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::build()
{
    base::build();
}

template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::clear()
{
    base::clear();
}

template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::flush()
{
    base::flush();
}

template <class GI, class FI, class SI, class PI, class MI> inline void
multipedigree<GI,FI,SI,PI,MI>::freeze()
{
    base::freeze();
}

//----------------------------------------------------------------------------
//
template <class GI, class FI, class SI, class PI, class MI> pedigree_id
multipedigree<GI,FI,SI,PI,MI>::allocate_pedigree(const string& name)
{
    return new pedigree_type(name, *this);
}

template <class GI, class FI, class SI, class PI, class MI> void
multipedigree<GI,FI,SI,PI,MI>::deallocate_pedigree(pedigree_id P)
{
    delete reinterpret_cast<pedigree_type*>(P);
}

} // End namespace MPED
} // End namespace SAGE
