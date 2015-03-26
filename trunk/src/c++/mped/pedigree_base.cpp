//============================================================================
//  File:       spbase1.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <mpfwd.h>
    #pragma hdrstop
#endif

#include "mped/spbase.h"

namespace SAGE {
namespace MPED {

static  uint    next_sp_index=0;
static  string  unconn("unconnected");

//----------------------------------------------------------------------------
//  Function:   ifind_family()
//
//  Purpose:    This member function returns an iterator from the set of
//              nuclear families, given the names of both parents.
//----------------------------------------------------------------------------
//
fset_iterator
pedigree_base::ifind_family(const string& parent1, const string& parent2) const
{
    member_base     m1(parent1, no_info());
    member_base     m2(parent2, no_info());
    family_base     f1(&m1, &m2, no_info());

    return my_families.find(&f1);
}


//----------------------------------------------------------------------------
//  Function:   ifind_family()
//
//  Purpose:    This member function returns a family id, given the member
//              id of both parents.
//----------------------------------------------------------------------------
//
fset_iterator
pedigree_base::ifind_family(member_id par1, member_id par2) const
{
    if (par1  &&  par2)
    {
        family_base     f1(par1, par2, no_info());
        
        return my_families.find(&f1);
    }
    else
    {
        return my_families.end();;
    }
}


//----------------------------------------------------------------------------
//  Function:   ifind_member()
//
//  Purpose:    This member function returns an iterator from the set of
//              members, given the name of that member.
//----------------------------------------------------------------------------
//
mset_iterator
pedigree_base::ifind_member(const string& name) const
{
    member_base     m1(name, no_info());

    return my_members.find(&m1);
}


//----------------------------------------------------------------------------
//  Function:   ifind_subped()
//
//  Purpose:    This member function returns an iterator from the list of
//              subpedigrees, given the name of that subpedigree.
//----------------------------------------------------------------------------
//
sset_iterator
pedigree_base::ifind_subped(const string& name) const
{
    subpedigree_base    s1(name, no_info());

    return my_subpeds.find(&s1);
}


//----------------------------------------------------------------------------
//  Function:   lookup_family()
//
//  Purpose:    This member function returns a family id, given the names of
//              both parents.
//----------------------------------------------------------------------------
//
family_id
pedigree_base::lookup_family(const string& parent1, const string& parent2) const
{
    member_base     m1(parent1, no_info());
    member_base     m2(parent2, no_info());
    family_base     f1(&m1, &m2, no_info());
    fset_iterator   it = my_families.find(&f1);

    return (it == my_families.end())  ?  family_id(0)  :  (*it);
}


//----------------------------------------------------------------------------
//  Function:   lookup_family()
//
//  Purpose:    This member function returns a family id, given the member
//              id of both parents.
//----------------------------------------------------------------------------
//
family_id
pedigree_base::lookup_family(member_id par1, member_id par2) const
{
    if (par1  &&  par2)
    {
        family_base     f1(par1, par2, no_info());
        fset_iterator   it = my_families.find(&f1);

        return (it == my_families.end())  ?  family_id(0)  :  (*it);
    }
    else
    {
        return 0;
    }
}


//----------------------------------------------------------------------------
//  Function:   lookup_member()
//
//  Purpose:    This member function returns a member id, given the name of
//              a member.
//----------------------------------------------------------------------------
//
member_id
pedigree_base::lookup_member(const string& name) const
{
    member_base     m1(name, no_info());
    mset_iterator   it = my_members.find(&m1);

    return (it == my_members.end())  ?  member_id(0)  :  (*it);
}


//----------------------------------------------------------------------------
//  Function:   lookup_subped()
//
//  Purpose:    This member function returns a subpedigree id, given the
//              name of that subpedigree.
//----------------------------------------------------------------------------
//
subped_id
pedigree_base::lookup_subped(const string& name) const
{
    subpedigree_base    s1(name, no_info());
    sset_iterator       it = my_subpeds.find(&s1);

    return (it == my_subpeds.end())  ?  subped_id(0)  :  (*it);
}

//----------------------------------------------------------------------------
//  Function:   
//
//  Purpose:    
//----------------------------------------------------------------------------
//
subped_id
pedigree_base::allocate_subped(const string& name, pedigree_id pped)
{
    return new subpedigree_base(name, pped);
}


family_id
pedigree_base::allocate_family()
{
    return new family_base();
}


member_id
pedigree_base::allocate_member(const string& name, SexCode x, const void*)
{
    return new member_base(name, x);
}


void
pedigree_base::deallocate_subped(subped_id S)
{
    delete S;
}


void
pedigree_base::deallocate_family(family_id F)
{
    delete F;
}


void
pedigree_base::deallocate_member(member_id M)
{
    delete M;
}

//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_begin(const string& m) const
{
    member_id   M = lookup_member(m);

    if (M  &&  M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.first);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given an identifier of the member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_begin(const member_base& m) const
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.first);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given a cursor pointing to the 
//              member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_begin(const member_const_iterator& m) const
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.first);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_end(const string& m) const
{
    member_id   M = lookup_member(m);

    if (M  &&  M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.second);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the
//              past-the-end mate of a marriage chain.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_end(const member_base& m) const
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.second);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the
//              past-the-end mate of a marriage chain.
//----------------------------------------------------------------------------
//
pedigree_base::mate_const_iterator
pedigree_base::mate_end(const member_const_iterator& m) const
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(M);

        return mate_const_iterator(range.second);
    }
    else
    {
        return mate_const_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the id of the relevant
//              nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_begin(const family_base& f) const
{
    family_id   F = const_cast<family_id>(&f);

    if (F->pedigree() == this)
    {
        return offspring_const_iterator(F->offspring());
    }
    else
    {
        return offspring_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given a cursor pointing to the
//              relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_begin(const family_const_iterator& f) const
{
    family_id   F = const_cast<family_id>(&(*f));

    if (F->pedigree() == this)
    {
        return offspring_const_iterator(F->offspring());
    }
    else
    {
        return offspring_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the names of both parents
//              of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_begin(const string& p1, const string& p2) const
{
    family_id   F = lookup_family(p1, p2);

    if (F)
    {
        return offspring_const_iterator(F->offspring());
    }
    else
    {
        return offspring_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the identifiers of both 
//              parents of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_begin
(const member_base& p1, const member_base& p2) const
{
    family_id   F = lookup_family(p1.name(), p2.name());

    if (F != 0)
    {
        return offspring_const_iterator(F->offspring());
    }
    else
    {
        return offspring_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given cursors pointing to both 
//              parents of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_begin
(const member_const_iterator& p1, const member_const_iterator& p2) const
{
    family_id   F = lookup_family(p1->name(), p2->name());

    if (F != 0)
    {
        return offspring_const_iterator(F->offspring());
    }
    else
    {
        return offspring_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_const_iterator
pedigree_base::offspring_end() const
{
    return offspring_const_iterator(0);
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a member, given that member's name.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_begin(const string& name) const
{
    member_id   M = lookup_member(name);
    family_id   F = (M) ? M->family() : 0;

    return parent_const_iterator(F);
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_begin(const member_base& m) const
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        return parent_const_iterator(M->family());
    }
    else
    {
        return parent_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_begin(const member_const_iterator& m) const
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        return parent_const_iterator(M->family());
    }
    else
    {
        return parent_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_begin(const family_base& f) const
{
    family_id   F = const_cast<family_id>(&f);

    if (F->pedigree() == this)
    {
        return parent_const_iterator(F);
    }
    else
    {
        return parent_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_begin(const family_const_iterator& f) const
{
    family_id   F = const_cast<family_id>(&(*f));

    if (F->pedigree() == this)
    {
        return parent_const_iterator(F);
    }
    else
    {
        return parent_const_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              parent sequence, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_const_iterator
pedigree_base::parent_end() const
{
    return parent_const_iterator(0);
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_begin(const string& M) const
{
    member_id   m = lookup_member(M);

    if (m  &&  m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.first, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_begin(const member_base& M) const
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.first, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_begin(const member_const_iterator& M) const
{
    member_id   m = const_cast<member_id>(&(*M));

    if (m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.first, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_end(const string& M) const
{
    member_id   m = lookup_member(M);

    if (m  &&  m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.second, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_end(const member_base& M) const
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.second, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_const_iterator
pedigree_base::progeny_end(const member_const_iterator& M) const
{
    member_id   m = const_cast<member_id>(&(*M));

    if (m->pedigree() == this)
    {
        const_mate_range  range = my_matechains.equal_range(m);

        return progeny_const_iterator(range.second, range.second);
    }
    else
    {
        return progeny_const_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given the name of one of the 
//              siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_const_iterator
pedigree_base::sibling_begin(const string& name) const
{
    if (name.size() > 0)
    {
        member_id   m1 = lookup_member(name);
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_const_iterator(m1, m2);
    }
    else
    {
        return sibling_const_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given a cursor pointing to one of
//              the siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_const_iterator
pedigree_base::sibling_begin(const member_base& M) const
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        member_id   m1 = m;
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_const_iterator(m1, m2);
    }
    else
    {
        return sibling_const_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given a cursor pointing to one of
//              the siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_const_iterator
pedigree_base::sibling_begin(const member_const_iterator& m) const
{
    if (m->pedigree() == this)
    {
        member_id   m1 = const_cast<member_id>(&(*m));
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_const_iterator(m1, m2);
    }
    else
    {
        return sibling_const_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              sibling chain.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_const_iterator
pedigree_base::sibling_end() const
{
    return sibling_const_iterator();
}

//----------------------------------------------------------------------------
//  Function:   index_swap_*()
//
//  Purpose:    These functions swap entries in the relevant indexing vectors.
//----------------------------------------------------------------------------
//
void
pedigree_base::subpedigree_index_swap(uint i, uint j)
{
    std::swap(my_subped_index[i], my_subped_index[j]);
    my_subped_index[i]->set_index(i);
    my_subped_index[j]->set_index(j);
}


void
pedigree_base::member_index_swap(uint i, uint j)
{
    std::swap(my_member_index[i], my_member_index[j]);
    my_member_index[i]->set_index(i);
    my_member_index[j]->set_index(j);
}


void
pedigree_base::family_index_swap(uint i, uint j)
{
    std::swap(my_family_index[i], my_family_index[j]);
    my_family_index[i]->set_index(i);
    my_family_index[j]->set_index(j);
}


//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_begin(const string& m)
{
    member_id   M = lookup_member(m);

    if (M  &&  M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.first);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given an identifier of the member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_begin(const member_base& m)
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.first);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given a cursor pointing to the 
//              member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_begin(const member_const_iterator& m)
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.first);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              mate of a marriage chain, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_end(const string& m)
{
    member_id   M = lookup_member(m);

    if (M  &&  M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.second);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the
//              past-the-end mate of a marriage chain.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_end(const member_base& m)
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.second);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   mate_end()
//
//  Purpose:    This member function returns a cursor pointing to the
//              past-the-end mate of a marriage chain.
//----------------------------------------------------------------------------
//
pedigree_base::mate_iterator
pedigree_base::mate_end(const member_const_iterator& m)
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(M);

        return mate_iterator(range.second);
    }
    else
    {
        return mate_iterator(my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the id of the relevant
//              nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_begin(const family_base& f)
{
    family_id   F = const_cast<family_id>(&f);

    if (F->pedigree() == this)
    {
        return offspring_iterator(F->offspring());
    }
    else
    {
        return offspring_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given a cursor pointing to the
//              relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_begin(const family_const_iterator& f)
{
    family_id   F = const_cast<family_id>(&(*f));

    if (F->pedigree() == this)
    {
        return offspring_iterator(F->offspring());
    }
    else
    {
        return offspring_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the names of both parents
//              of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_begin(const string& p1, const string& p2)
{
    family_id   F = lookup_family(p1, p2);

    if (F)
    {
        return offspring_iterator(F->offspring());
    }
    else
    {
        return offspring_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given the identifiers of both 
//              parents of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_begin
(const member_base& p1, const member_base& p2)
{
    family_id   F = lookup_family(p1.name(), p2.name());

    if (F != 0)
    {
        return offspring_iterator(F->offspring());
    }
    else
    {
        return offspring_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of an offspring chain, given cursors pointing to both 
//              parents of the relevant nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_begin
(const member_const_iterator& p1, const member_const_iterator& p2)
{
    family_id   F = lookup_family(p1->name(), p2->name());

    if (F != 0)
    {
        return offspring_iterator(F->offspring());
    }
    else
    {
        return offspring_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   offspring_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::offspring_iterator
pedigree_base::offspring_end()
{
    return offspring_iterator(0);
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a member, given that member's name.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_begin(const string& name)
{
    member_id   M = lookup_member(name);
    family_id   F = (M) ? M->family() : 0;

    return parent_iterator(F);
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_begin(const member_base& m)
{
    member_id   M = const_cast<member_id>(&m);

    if (M->pedigree() == this)
    {
        return parent_iterator(M->family());
    }
    else
    {
        return parent_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_begin(const member_const_iterator& m)
{
    member_id   M = const_cast<member_id>(&(*m));

    if (M->pedigree() == this)
    {
        return parent_iterator(M->family());
    }
    else
    {
        return parent_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_begin(const family_base& f)
{
    family_id   F = const_cast<family_id>(&f);

    if (F->pedigree() == this)
    {
        return parent_iterator(F);
    }
    else
    {
        return parent_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              parent of a given nuclear family.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_begin(const family_const_iterator& f)
{
    family_id   F = const_cast<family_id>(&(*f));

    if (F->pedigree() == this)
    {
        return parent_iterator(F);
    }
    else
    {
        return parent_iterator(0);
    }
}


//----------------------------------------------------------------------------
//  Function:   parent_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              parent sequence, given the name of a member.
//----------------------------------------------------------------------------
//
pedigree_base::parent_iterator
pedigree_base::parent_end()
{
    return parent_iterator(0);
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_begin(const string& M)
{
    member_id   m = lookup_member(M);

    if (m  &&  m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.first, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_begin(const member_base& M)
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.first, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_begin()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_begin(const member_const_iterator& M)
{
    member_id   m = const_cast<member_id>(&(*M));

    if (m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.first, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_end(const string& M)
{
    member_id   m = lookup_member(M);

    if (m  &&  m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.second, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_end(const member_base& M)
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.second, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   progeny_end()
//
//  Purpose:    This member function returns a past-the-end cursor for a
//              total offspring chain.
//----------------------------------------------------------------------------
//
pedigree_base::progeny_iterator
pedigree_base::progeny_end(const member_const_iterator& M)
{
    member_id   m = const_cast<member_id>(&(*M));

    if (m->pedigree() == this)
    {
        mate_range  range = my_matechains.equal_range(m);

        return progeny_iterator(range.second, range.second);
    }
    else
    {
        return progeny_iterator(my_matechains.end(), my_matechains.end());
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given the name of one of the 
//              siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_iterator
pedigree_base::sibling_begin(const string& name)
{
    if (name.size() > 0)
    {
        member_id   m1 = lookup_member(name);
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_iterator(m1, m2);
    }
    else
    {
        return sibling_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given a cursor pointing to one of
//              the siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_iterator
pedigree_base::sibling_begin(const member_base& M)
{
    member_id   m = const_cast<member_id>(&M);

    if (m->pedigree() == this)
    {
        member_id   m1 = m;
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_iterator(m1, m2);
    }
    else
    {
        return sibling_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_begin()
//
//  Purpose:    This member function returns a cursor pointing to the first
//              member of a sibling chain, given a cursor pointing to one of
//              the siblings.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_iterator
pedigree_base::sibling_begin(const member_const_iterator& m)
{
    if (m->pedigree() == this)
    {
        member_id   m1 = const_cast<member_id>(&(*m));
        family_id   f1 = (m1)  ?  m1->family()     :  0;
        member_id   m2 = (f1)  ?  f1->offspring()  :  0;

        return sibling_iterator(m1, m2);
    }
    else
    {
        return sibling_iterator();
    }
}


//----------------------------------------------------------------------------
//  Function:   sibling_end()
//
//  Purpose:    This member function returns a past-the-end cursor for an
//              sibling chain.
//----------------------------------------------------------------------------
//
pedigree_base::sibling_iterator
pedigree_base::sibling_end()
{
    return sibling_iterator();
}

//============================================================================
//  IMPLEMENTATION: pedigree_base
//============================================================================
//
//- Canonical members of "pedigree_base".
//
pedigree_base::pedigree_base(const string& name)
  : my_fcount(0), my_mcount(0), my_scount(0), my_ucount(0), 
    my_name(name), 
    my_readiness(false), my_fstate(false),
    my_unconnecteds(0), my_mped(0)
{
    my_unconnecteds = pedigree_base::allocate_subped(unconn, this);
}


pedigree_base::pedigree_base(const string& name, multipedigree_pointer mp)
  : my_fcount(0), my_mcount(0), my_scount(0), my_ucount(0), 
    my_name(name), 
    my_readiness(false), my_fstate(false),
    my_unconnecteds(0), my_mped(mp)
{
    my_unconnecteds = pedigree_base::allocate_subped(unconn, this);
}


pedigree_base::~pedigree_base()
{
    if (my_unconnecteds != 0)
    {
        pedigree_base::deallocate_subped(my_unconnecteds);
    }

    if (my_members.size() > 0)
    {
        clear();
    }
}


//----------------------------------------------------------------------------
//  Function:   add_lineage()
//
//  Purpose:    This member function creates a lineage meta-relationship
//              given the name of the child and the name of one parent.
//              It performs validation of all names.
//----------------------------------------------------------------------------
//
void
pedigree_base::add_lineage(const string& child, const string& parent)
{
    if (my_fstate == false)
    {
        my_readiness = !my_builder.add_lineage(child, parent);
    }
}


//----------------------------------------------------------------------------
//  Function:   add_lineage()
//
//  Purpose:    This member function creates a lineage meta-relationship
//              given the name of the child and the names of both parents.
//              It performs validation of all names.
//----------------------------------------------------------------------------
//
void
pedigree_base::add_lineage
(const string& child, const string& parent1, const string& parent2)
{
    if (my_fstate == false)
    {
        my_readiness = !my_builder.add_lineage(child, parent1, parent2);
    }
}


//----------------------------------------------------------------------------
//  Function:   add_marriage()
//
//  Purpose:    This member function adds a marriage meta-relationship, given
//              the names of both mates.
//----------------------------------------------------------------------------
//
void
pedigree_base::add_marriage(const string& mate1, const string& mate2)
{
    if (my_fstate == false)
    {
        my_readiness = !my_builder.add_marriage(mate1, mate2);
    }
}


//----------------------------------------------------------------------------
//  Function:   add_member()
//
//  Purpose:    This member function adds a new member to the pedigree.
//----------------------------------------------------------------------------
//
void
pedigree_base::add_member(const string& name, SexCode s)
{
    if (my_fstate == false)
    {
        add_member(name, s, 0);
    }
}


//----------------------------------------------------------------------------
//  Function:   add_sibship()
//
//  Purpose:    This member function adds a sibling meta-relationship, given
//              the names of both siblings.
//----------------------------------------------------------------------------
//
void
pedigree_base::add_sibship(const string& sib1, const string& sib2)
{
    if (my_fstate == false)
    {
        my_readiness = !my_builder.add_sibship(sib1, sib2);
    }
}


//----------------------------------------------------------------------------
//  Function:   build()
//
//  Purpose:    This member function constructs as much of a pedigree as 
//              possible given the meta-relationship information the pedigree
//              object currently possesses.
//----------------------------------------------------------------------------
//
void
pedigree_base::build()
{
    if (my_fstate == false)
    {
        my_builder.build_pedigree(*this);
        
        my_readiness = true;
    }
}


//----------------------------------------------------------------------------
//  Function:   clear()
//
//  Purpose:    This member function deletes all meta-relationships, members, 
//              families, and subpedigrees; it also clears all STL containers
//              used to hold and/or represent this information.
//----------------------------------------------------------------------------
//
void
pedigree_base::clear()
{
    //- Clear all remaining meta-relationship information.
    //
    flush_build_info();

    //- Delete the family object associated with each family ID in the
    //  set of families.
    //
    family_set::iterator    fsf = my_families.begin();
    family_set::iterator    fsl = my_families.end();

    for (;  fsf != fsl;  ++fsf)
    {
        deallocate_family(*fsf);
    }

    //- Delete the member object associated with each member ID in the
    //  set of members.
    //
    member_set::iterator    msf = my_members.begin();
    member_set::iterator    msl = my_members.end();

    for (;  msf != msl;  ++msf)
    {
        deallocate_member(*msf);
    }

    //- Delete the member object associated with each subped ID in the
    //  set of subpedigrees.
    //
    subped_set::iterator   slf = my_subpeds.begin();
    subped_set::iterator   sll = my_subpeds.end();

    for (;  slf != sll;  ++slf)
    {
        deallocate_subped(*slf);
    }

    //- Clear the sets/lists.
    //
    my_families.clear();
    my_members.clear();
    my_subpeds.clear();

    //- Allow new meta-relationship entries.
    //
    my_fstate = false;
}


//----------------------------------------------------------------------------
//  Function:   flush_build_info()
//
//  Purpose:    This member function deletes all meta-relationship information
//              and clears all STL containers used to hold and/or represent 
//              this information.
//----------------------------------------------------------------------------
//
void
pedigree_base::flush_build_info()
{
    my_builder.flush_build_info();
}


//----------------------------------------------------------------------------
//  Function:   freeze()
//
//  Purpose:    This member function freezes the state of the pedigree so
//              that new meta-relationship information cannot be added.
//----------------------------------------------------------------------------
//
void
pedigree_base::freeze()
{
    my_fstate = true;
}


//----------------------------------------------------------------------------
//  Function:   init()
//
//  Purpose:    This function allocates a subpedigree object for 
//              unconnected members.  This function should only be called
//              by constructors of the derived interface class.
//----------------------------------------------------------------------------
//
void
pedigree_base::init()
{
    string  uname;

    uname.assign(my_name).append(":unconnected");
    my_unconnecteds = allocate_subped(uname, this);
}


//----------------------------------------------------------------------------
//  Function:   clear_all()
//
//  Purpose:    This function deallocates the unconnected subpedigree object 
//              and clears the pedigree.
//----------------------------------------------------------------------------
//
void
pedigree_base::clear_all()
{
    pedigree_base::deallocate_subped(my_unconnecteds);
    my_unconnecteds = 0;
    clear();
}


//----------------------------------------------------------------------------
//  Function:   add_family()
//
//  Purpose:    This member function creates a new nuclear family and adds
//              it to the set of nuclear families belonging to this pedigree.
//----------------------------------------------------------------------------
//
family_id
pedigree_base::add_family(member_id P1, member_id P2, member_id C)
{
    family_id   id = allocate_family();

    id->reset_links(my_unconnecteds, P1, P2, C);
    my_families.insert(id);
    my_matechains.insert( mate_pair(P1, mate_info_base(P2, id)) );
    my_matechains.insert( mate_pair(P2, mate_info_base(P1, id)) );
    my_readiness = false;

    return id;
}


//----------------------------------------------------------------------------
//  Function:   add_member()
//
//  Purpose:    This member function creates a new member and adds
//              it to the set of members belonging to this pedigree.
//----------------------------------------------------------------------------
//
member_id
pedigree_base::add_member(const string& name, SexCode gender, const void* pinfo)
{
    member_id   id = member_find(name);

    if(id)
    {
      // Member has already been added
      if(gender != id->my_gender && !is_sex_unknown(gender))
      {
        // Do not down-grade known-sex to assumed-sex
        if ((id->is_male()   && gender == SEX_XMALE)  ||
            (id->is_female() && gender == SEX_XFEMALE)   )
        {
            gender = id->my_gender;
        }

        // Warn if member is trans-gendered
        if (id->is_male() && is_female(gender))
        {
            my_builder.my_errors.push_back(error_info(error_info::bad_gender, name, "male", "female"));
        }

        else if (id->is_female() && is_male(gender))
        {
            my_builder.my_errors.push_back(error_info(error_info::bad_gender, name, "female", "male"));
        }

        // Assign new gender
        id->set_sex(gender);
      }
    }
    else
    {
        id = allocate_member(name, gender, pinfo);
        id->set_subped(my_unconnecteds);
        my_members.insert(id);
        my_readiness = false;
    }

    return id;
}


//----------------------------------------------------------------------------
//  Function:   add_subped()
//
//  Purpose:    This member function creates a new member and adds
//              it to the set of members belonging to this pedigree.
//----------------------------------------------------------------------------
//
subped_id
pedigree_base::add_subped()
{
    char        cname[16];
    string      name;
    subped_id   id;

    sprintf(cname, ":%d", ++next_sp_index);
    name.assign(my_name).append(cname);
    id = allocate_subped(name, this);
    my_subpeds.insert(id);
    my_readiness = false;

    return id;
}

} // End namespace MPED
} // End namespace SAGE
