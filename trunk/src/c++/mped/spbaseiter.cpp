//============================================================================
//  File:       spbaseiter.cpp
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

//----------------------------------------------------------------------------
//  Class:      cursor_proxy
//
//  Purpose:    Provides access to pedigree_base functions without requiring
//              large numbers of friendships.
//----------------------------------------------------------------------------
//
class cursor_proxy
{
  public:
    static  family_id       origin(member_id M);
    static  member_id       parent1(member_id M);
    static  member_id       parent2(member_id M);
    static  pedigree_id     pedigree(member_id M);
    static  member_id       siblings(member_id M);
    static  subped_id       subpedigree(member_id M);

    static  member_id       offspring(family_id F);
    static  member_id       parent1(family_id F);
    static  member_id       parent2(family_id F);
    static  pedigree_id     pedigree(family_id F);
    static  subped_id       subpedigree(family_id F);

    static  pedigree_id     pedigree(subped_id S);

    static  family_base_iterator   family_begin(pedigree_id P, subped_id S);
    static  family_base_iterator   family_end(pedigree_id P);
    static  family_base_iterator   family_begin(subped_id P, subped_id S);
    static  family_base_iterator   family_end(subped_id P);
    static  member_base_iterator   member_begin(pedigree_id P, subped_id S);
    static  member_base_iterator   member_end(pedigree_id P);
    static  member_base_iterator   member_begin(subped_id P, subped_id S);
    static  member_base_iterator   member_end(subped_id P);

    static  const_mmap_iterator     mate_end(member_id);
};


//============================================================================
//  IMPLEMENTATION: cursor_proxy
//============================================================================
//
inline family_id
cursor_proxy::origin(member_id M)
{
    return M->family();
}

inline member_id
cursor_proxy::parent1(member_id M)
{
    return (M->family())  ?  M->family()->parent1()  :  member_id(0);
}

inline member_id
cursor_proxy::parent2(member_id M)
{
    return (M->family())  ?  M->family()->parent2()  :  member_id(0);
}

inline pedigree_id
cursor_proxy::pedigree(member_id M)
{
    return M->pedigree();
}

inline member_id
cursor_proxy::siblings(member_id M)
{
    return M->siblings();
}

inline subped_id
cursor_proxy::subpedigree(member_id M)
{
    return M->subpedigree();
}

//----------
//
inline member_id
cursor_proxy::offspring(family_id F)
{
    return F->offspring();
}

inline member_id
cursor_proxy::parent1(family_id F)
{
    return F->parent1();
}

inline member_id
cursor_proxy::parent2(family_id F)
{
    return F->parent2();
}

inline pedigree_id
cursor_proxy::pedigree(family_id F)
{
    return F->pedigree();
}

inline subped_id
cursor_proxy::subpedigree(family_id F)
{
    return F->subpedigree();
}

//----------
//
inline pedigree_id
cursor_proxy::pedigree(subped_id S)
{
    return S->pedigree();
}

//----------
//
inline family_base_iterator
cursor_proxy::family_begin(pedigree_id P, subped_id S)
{
    return P->family_begin(*S);
}

inline family_base_iterator
cursor_proxy::family_end(pedigree_id P)
{
    return P->family_end();
}

inline family_base_iterator
cursor_proxy::family_begin(subped_id P, subped_id S)
{
    return P->pedigree()->family_begin(*S);
}

inline family_base_iterator
cursor_proxy::family_end(subped_id P)
{
    return P->pedigree()->family_end();
}

//----------
//
inline member_base_iterator
cursor_proxy::member_begin(pedigree_id P, subped_id S)
{
    return P->member_begin(*S);
}

inline member_base_iterator
cursor_proxy::member_end(pedigree_id P)
{
    return P->member_end();
}

inline member_base_iterator
cursor_proxy::member_begin(subped_id P, subped_id S)
{
    return P->pedigree()->member_begin(*S);
}

inline member_base_iterator
cursor_proxy::member_end(subped_id P)
{
    return P->pedigree()->member_end();
}

inline const_mmap_iterator
cursor_proxy::mate_end(member_id M)
{
    return M->pedigree()->my_matechains.end();
}


//============================================================================
//  IMPLEMENTATION: offspring_base_iterator
//============================================================================
//
void
offspring_base_iterator::advance()
{
    if (my_data)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
}


//============================================================================
//  IMPLEMENTATION: offspring_base_const_iterator
//============================================================================
//
void
offspring_base_const_iterator::advance()
{
    if (my_data)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
}


//============================================================================
//  IMPLEMENTATION: parent_base_iterator
//============================================================================
//
void
parent_base_iterator::init()
{
    if (my_family)
    {
        my_data = cursor_proxy::parent1(my_family);
    }
}


void
parent_base_iterator::advance()
{
    if (my_family  &&  my_data)
    {
        if (my_data == cursor_proxy::parent1(my_family))
        {
            my_data = cursor_proxy::parent2(my_family);
        }
        else if (my_data == cursor_proxy::parent2(my_family))
        {
            my_data   = 0;
            my_family = 0;
        }
    }
}


//============================================================================
//  IMPLEMENTATION: parent_base_const_iterator
//============================================================================
//
void
parent_base_const_iterator::init()
{
    if (my_family)
    {
        my_data = cursor_proxy::parent1(my_family);
    }
}


void
parent_base_const_iterator::advance()
{
    if (my_family  &&  my_data)
    {
        if (my_data == cursor_proxy::parent1(my_family))
        {
            my_data = cursor_proxy::parent2(my_family);
        }
        else if (my_data == cursor_proxy::parent2(my_family))
        {
            my_data   = 0;
            my_family = 0;
        }
    }
}


//============================================================================
//  IMPLEMENTATION: progeny_base_iterator
//============================================================================
//
void
progeny_base_iterator::init()
{
    for (;  my_iter != my_limit;  ++my_iter)
    {
        my_data = cursor_proxy::offspring( &(my_iter->second.family()) );

        if (my_data != 0)
        {
            break;
        }
    }
}


void
progeny_base_iterator::advance()
{
    if (my_iter != my_limit)
    {
        my_data = cursor_proxy::siblings(my_data);

        while (my_data == 0)
        {
            ++my_iter;

            if (my_iter != my_limit)
            {
                my_data = cursor_proxy::offspring( &(my_iter->second.family()) );
            }
            else
            {
                break;
            }
        }
    }
}


//============================================================================
//  IMPLEMENTATION: progeny_base_const_iterator
//============================================================================
//
void
progeny_base_const_iterator::init()
{
    for (;  my_iter != my_limit;  ++my_iter)
    {
        family_id   ftmp;
        
        ftmp    = const_cast<family_id>( &(my_iter->second.family()) );
        my_data = cursor_proxy::offspring(ftmp);

        if (my_data != 0)
        {
            break;
        }
    }
}


void
progeny_base_const_iterator::advance()
{
    if (my_iter != my_limit)
    {
        my_data = cursor_proxy::siblings(my_data);

        while (my_data == 0)
        {
            ++my_iter;

            if (my_iter != my_limit)
            {
                family_id   ftmp;
                
                ftmp    = const_cast<family_id>( &(my_iter->second.family() ));
                my_data = cursor_proxy::offspring(ftmp);
            }
            else
            {
                break;
            }
        }
    }
}


//============================================================================
//  IMPLEMENTATION: sibling_base_iterator
//============================================================================
//
void
sibling_base_iterator::advance()
{
    if (my_data)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
    if (my_data  &&  my_data == my_self)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
}


//============================================================================
//  IMPLEMENTATION: sibling_base_const_iterator
//============================================================================
//
void
sibling_base_const_iterator::advance()
{
    if (my_data)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
    if (my_data  &&  my_data == my_self)
    {
        my_data = cursor_proxy::siblings(my_data);
    }
}

} // End namespace MPED
} // End namespace SAGE
