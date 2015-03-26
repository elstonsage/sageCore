#ifndef _MPITERATOR_HPP
#define _MPITERATOR_HPP

//============================================================================
//  File:       mpiterator.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.90
//
//  Notes:      
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mped/mpbase.h"

SAGE_NS_BEGIN

//----------------------------------------------------------------------------
//  Class:      pedigree_iterator<GI,FI,SI,PI,MI>
//
//  Purpose:    This class performs iteration over the pedigrees in
//              a multipedigree.
//----------------------------------------------------------------------------
//
template <class GI, class FI, class SI, class PI, class MI>
class pedigree_iterator
{
  public:
    friend  class multipedigree<GI,FI,SI,PI,MI>;

    typedef pedigree_iterator<GI,FI,SI,PI,MI>   this_type;
    typedef pedigree<GI,FI,SI,PI,MI>&           reference;
    typedef pedigree<GI,FI,SI,PI,MI>*           pointer;
    typedef pedigree_map::iterator              cursor_type;
    typedef const pedigree_map::iterator&       cursor_reference;
        

  public:
    pedigree_iterator();

    bool            operator ==(const pedigree_iterator& i) const;
    bool            operator !=(const pedigree_iterator& i) const;
    pointer         operator ->() const;
    reference       operator *() const;

    this_type&      operator ++();
    this_type       operator ++(int);

    cursor_reference    cursor() const;

  private:
    cursor_type     my_cursor;

    pedigree_iterator(cursor_reference c);
};


//============================================================================
//  IMPLEMENTATION: pedigree_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_iterator<GI,FI,SI,PI,MI>::pedigree_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_iterator<GI,FI,SI,PI,MI>::pedigree_iterator(cursor_reference c)
  : my_cursor(c)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_iterator<GI,FI,SI,PI,MI>::operator ==(const pedigree_iterator& i) const
{
    return my_cursor == i.my_cursor;
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_iterator<GI,FI,SI,PI,MI>::operator !=(const pedigree_iterator& i) const
{
    return my_cursor != i.my_cursor;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>::pointer
pedigree_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return static_cast<pointer>( my_cursor->second.data() );
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>::reference
pedigree_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return *static_cast<pointer>( my_cursor->second.data() );
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>&
pedigree_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    ++my_cursor;
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>
pedigree_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    cursor_type     tmp(*this);

    ++my_cursor
    return tmp;
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>::cursor_reference
pedigree_iterator<GI,FI,SI,PI,MI>::cursor() const
{
    return my_cursor;
}

SAGE_NS_END
#endif  //- _MPITERATOR_HPP
