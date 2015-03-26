//============================================================================
//  File:       spiter.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//  IMPLEMENTATION: family_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
namespace SAGE {
namespace MPED {

template <class GI, class FI, class SI, class PI, class MI> inline
family_iterator<GI,FI,SI,PI,MI>::family_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
family_iterator<GI,FI,SI,PI,MI>::family_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
family_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
family_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
family_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_iterator<GI,FI,SI,PI,MI>::pointer
family_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_iterator<GI,FI,SI,PI,MI>::reference
family_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_iterator<GI,FI,SI,PI,MI>::reference
family_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>&
family_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>&
family_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>&
family_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>&
family_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>
family_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_iterator<GI,FI,SI,PI,MI>
family_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: family_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
family_const_iterator<GI,FI,SI,PI,MI>::family_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
family_const_iterator<GI,FI,SI,PI,MI>::family_const_iterator
(const family_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
family_const_iterator<GI,FI,SI,PI,MI>::family_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
family_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
family_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
family_const_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_const_iterator<GI,FI,SI,PI,MI>::pointer
family_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_const_iterator<GI,FI,SI,PI,MI>::reference
family_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename family_const_iterator<GI,FI,SI,PI,MI>::reference
family_const_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>&
family_const_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>&
family_const_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>&
family_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>&
family_const_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>
family_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
family_const_iterator<GI,FI,SI,PI,MI>
family_const_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: mate_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
mate_iterator<GI,FI,SI,PI,MI>::mate_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
mate_iterator<GI,FI,SI,PI,MI>::mate_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
mate_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
mate_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_iterator<GI,FI,SI,PI,MI>::pointer
mate_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_iterator<GI,FI,SI,PI,MI>::reference
mate_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
mate_iterator<GI,FI,SI,PI,MI>&
mate_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
mate_iterator<GI,FI,SI,PI,MI>
mate_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: mate_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
mate_const_iterator<GI,FI,SI,PI,MI>::mate_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
mate_const_iterator<GI,FI,SI,PI,MI>::mate_const_iterator
(const mate_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
mate_const_iterator<GI,FI,SI,PI,MI>::mate_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
mate_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
mate_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_const_iterator<GI,FI,SI,PI,MI>::pointer
mate_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_const_iterator<GI,FI,SI,PI,MI>::reference
mate_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
mate_const_iterator<GI,FI,SI,PI,MI>&
mate_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
mate_const_iterator<GI,FI,SI,PI,MI>
mate_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: member_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
member_iterator<GI,FI,SI,PI,MI>::member_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
member_iterator<GI,FI,SI,PI,MI>::member_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
member_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_iterator<GI,FI,SI,PI,MI>::pointer
member_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_iterator<GI,FI,SI,PI,MI>::reference
member_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_iterator<GI,FI,SI,PI,MI>::reference
member_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>&
member_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>&
member_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>&
member_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>&
member_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>
member_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_iterator<GI,FI,SI,PI,MI>
member_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: member_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
member_const_iterator<GI,FI,SI,PI,MI>::member_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
member_const_iterator<GI,FI,SI,PI,MI>::member_const_iterator
(const member_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
member_const_iterator<GI,FI,SI,PI,MI>::member_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
member_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
member_const_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_const_iterator<GI,FI,SI,PI,MI>::pointer
member_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_const_iterator<GI,FI,SI,PI,MI>::reference
member_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename member_const_iterator<GI,FI,SI,PI,MI>::reference
member_const_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>&
member_const_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>&
member_const_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>&
member_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>&
member_const_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>
member_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
member_const_iterator<GI,FI,SI,PI,MI>
member_const_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: offspring_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
offspring_iterator<GI,FI,SI,PI,MI>::offspring_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
offspring_iterator<GI,FI,SI,PI,MI>::offspring_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
offspring_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
offspring_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename offspring_iterator<GI,FI,SI,PI,MI>::pointer
offspring_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename offspring_iterator<GI,FI,SI,PI,MI>::reference
offspring_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
offspring_iterator<GI,FI,SI,PI,MI>&
offspring_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
offspring_iterator<GI,FI,SI,PI,MI>
offspring_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: offspring_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
offspring_const_iterator<GI,FI,SI,PI,MI>::offspring_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
offspring_const_iterator<GI,FI,SI,PI,MI>::offspring_const_iterator
(const offspring_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
offspring_const_iterator<GI,FI,SI,PI,MI>::offspring_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
offspring_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
offspring_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename offspring_const_iterator<GI,FI,SI,PI,MI>::pointer
offspring_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename offspring_const_iterator<GI,FI,SI,PI,MI>::reference
offspring_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
offspring_const_iterator<GI,FI,SI,PI,MI>&
offspring_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
offspring_const_iterator<GI,FI,SI,PI,MI>
offspring_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: parent_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
parent_iterator<GI,FI,SI,PI,MI>::parent_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
parent_iterator<GI,FI,SI,PI,MI>::parent_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
parent_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
parent_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename parent_iterator<GI,FI,SI,PI,MI>::pointer
parent_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename parent_iterator<GI,FI,SI,PI,MI>::reference
parent_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
parent_iterator<GI,FI,SI,PI,MI>&
parent_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
parent_iterator<GI,FI,SI,PI,MI>
parent_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

//============================================================================
//  IMPLEMENTATION: parent_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
parent_const_iterator<GI,FI,SI,PI,MI>::parent_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
parent_const_iterator<GI,FI,SI,PI,MI>::parent_const_iterator
(const parent_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
parent_const_iterator<GI,FI,SI,PI,MI>::parent_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
parent_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
parent_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename parent_const_iterator<GI,FI,SI,PI,MI>::pointer
parent_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename parent_const_iterator<GI,FI,SI,PI,MI>::reference
parent_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
parent_const_iterator<GI,FI,SI,PI,MI>&
parent_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
parent_const_iterator<GI,FI,SI,PI,MI>
parent_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: pedigree_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_iterator<GI,FI,SI,PI,MI>::pedigree_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_iterator<GI,FI,SI,PI,MI>::pedigree_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_iterator<GI,FI,SI,PI,MI>::pointer
pedigree_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_iterator<GI,FI,SI,PI,MI>::reference
pedigree_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_iterator<GI,FI,SI,PI,MI>::reference
pedigree_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>&
pedigree_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>&
pedigree_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>&
pedigree_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>&
pedigree_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>
pedigree_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_iterator<GI,FI,SI,PI,MI>
pedigree_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: pedigree_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_const_iterator<GI,FI,SI,PI,MI>::pedigree_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_const_iterator<GI,FI,SI,PI,MI>::pedigree_const_iterator
(const pedigree_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
pedigree_const_iterator<GI,FI,SI,PI,MI>::pedigree_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_const_iterator<GI,FI,SI,PI,MI>::pointer
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_const_iterator<GI,FI,SI,PI,MI>::reference
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename pedigree_const_iterator<GI,FI,SI,PI,MI>::reference
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>&
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>&
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>&
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>&
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
pedigree_const_iterator<GI,FI,SI,PI,MI>
pedigree_const_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: progeny_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
progeny_iterator<GI,FI,SI,PI,MI>::progeny_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
progeny_iterator<GI,FI,SI,PI,MI>::progeny_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
progeny_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
progeny_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename progeny_iterator<GI,FI,SI,PI,MI>::pointer
progeny_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename progeny_iterator<GI,FI,SI,PI,MI>::reference
progeny_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
progeny_iterator<GI,FI,SI,PI,MI>&
progeny_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
progeny_iterator<GI,FI,SI,PI,MI>
progeny_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: progeny_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
progeny_const_iterator<GI,FI,SI,PI,MI>::progeny_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
progeny_const_iterator<GI,FI,SI,PI,MI>::progeny_const_iterator
(const progeny_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
progeny_const_iterator<GI,FI,SI,PI,MI>::progeny_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
progeny_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
progeny_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename progeny_const_iterator<GI,FI,SI,PI,MI>::pointer
progeny_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename progeny_const_iterator<GI,FI,SI,PI,MI>::reference
progeny_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
progeny_const_iterator<GI,FI,SI,PI,MI>&
progeny_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
progeny_const_iterator<GI,FI,SI,PI,MI>
progeny_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: sibling_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
sibling_iterator<GI,FI,SI,PI,MI>::sibling_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
sibling_iterator<GI,FI,SI,PI,MI>::sibling_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
sibling_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
sibling_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename sibling_iterator<GI,FI,SI,PI,MI>::pointer
sibling_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename sibling_iterator<GI,FI,SI,PI,MI>::reference
sibling_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
sibling_iterator<GI,FI,SI,PI,MI>&
sibling_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
sibling_iterator<GI,FI,SI,PI,MI>
sibling_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: sibling_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
sibling_const_iterator<GI,FI,SI,PI,MI>::sibling_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
sibling_const_iterator<GI,FI,SI,PI,MI>::sibling_const_iterator
(const sibling_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
sibling_const_iterator<GI,FI,SI,PI,MI>::sibling_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
sibling_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
sibling_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename sibling_const_iterator<GI,FI,SI,PI,MI>::pointer
sibling_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename sibling_const_iterator<GI,FI,SI,PI,MI>::reference
sibling_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
sibling_const_iterator<GI,FI,SI,PI,MI>&
sibling_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
sibling_const_iterator<GI,FI,SI,PI,MI>
sibling_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: subpedigree_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree_iterator<GI,FI,SI,PI,MI>::subpedigree_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree_iterator<GI,FI,SI,PI,MI>::subpedigree_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_iterator<GI,FI,SI,PI,MI>::pointer
subpedigree_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_iterator<GI,FI,SI,PI,MI>::reference
subpedigree_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_iterator<GI,FI,SI,PI,MI>::reference
subpedigree_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>&
subpedigree_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>&
subpedigree_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>&
subpedigree_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>&
subpedigree_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>
subpedigree_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_iterator<GI,FI,SI,PI,MI>
subpedigree_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: subpedigree_const_iterator<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree_const_iterator<GI,FI,SI,PI,MI>::subpedigree_const_iterator()
{}

template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree_const_iterator<GI,FI,SI,PI,MI>::subpedigree_const_iterator
(const subpedigree_iterator<GI,FI,SI,PI,MI>& i)
  : base_type(i)
{}

template <class GI, class FI, class SI, class PI, class MI> inline
subpedigree_const_iterator<GI,FI,SI,PI,MI>::subpedigree_const_iterator
(const base_type& i)
  : base_type(i)
{}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator ==(const this_type& i) const
{
    return base_type::operator==(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator !=(const this_type& i) const
{
    return base_type::operator!=(i);
}

template <class GI, class FI, class SI, class PI, class MI> inline bool
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator <(const this_type& i) const
{
    return base_type::operator<(i);
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_const_iterator<GI,FI,SI,PI,MI>::pointer
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator ->() const
{
    return reinterpret_cast<pointer>(base_type::operator->()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_const_iterator<GI,FI,SI,PI,MI>::reference
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator *() const
{
    return reinterpret_cast<reference>(base_type::operator*()); 
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename subpedigree_const_iterator<GI,FI,SI,PI,MI>::reference
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator [](difference_type n) const
{
    return reinterpret_cast<reference>(base_type::operator[](n)); 
}

//----------
//
template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>&
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator +=(difference_type n)
{
    base_type::operator+=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>&
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator -=(difference_type n)
{
    base_type::operator-=(n);
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>&
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator ++()
{
    base_type::operator++();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>&
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator --()
{
    base_type::operator--();
    return *this;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator ++(int)
{
    this_type   tmp(*this);

    base_type::operator++((int)0);
    return tmp;
}

template <class GI, class FI, class SI, class PI, class MI> inline 
subpedigree_const_iterator<GI,FI,SI,PI,MI>
subpedigree_const_iterator<GI,FI,SI,PI,MI>::operator --(int)
{
    this_type   tmp(*this);

    base_type::operator--((int)0);
    return tmp;
}

} // End namespace MPED
} // End namespace SAGE
