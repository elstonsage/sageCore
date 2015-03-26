//============================================================================
//  File:       spbaseiter.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//
//============================================================================
//  IMPLEMENTATION: family_base_iterator
//============================================================================
//
namespace SAGE {
namespace MPED {

inline
family_base_iterator::family_base_iterator() 
  : my_iter()
{}

inline
family_base_iterator::family_base_iterator(const cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
family_base_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
family_base_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
family_base_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline family_base_iterator::pointer
family_base_iterator::operator ->() const
{
    return *my_iter;
}

inline family_base_iterator::reference
family_base_iterator::operator *() const
{
    return **my_iter;
}

inline family_base_iterator::reference
family_base_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline family_base_iterator&
family_base_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline family_base_iterator&
family_base_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline family_base_iterator&
family_base_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline family_base_iterator&
family_base_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline family_base_iterator
family_base_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline family_base_iterator
family_base_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: family_base_const_iterator
//============================================================================
//
inline
family_base_const_iterator::family_base_const_iterator() 
  : my_iter()
{}

inline
family_base_const_iterator::family_base_const_iterator
(const family_base_iterator& iter) 
  : my_iter(iter.my_iter)
{}

inline
family_base_const_iterator::family_base_const_iterator
(const cursor_type& iter) 
  : my_iter(iter)
{}

inline
family_base_const_iterator::family_base_const_iterator
(const const_cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
family_base_const_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
family_base_const_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
family_base_const_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline family_base_const_iterator::pointer
family_base_const_iterator::operator ->() const
{
    return *my_iter;
}

inline family_base_const_iterator::reference
family_base_const_iterator::operator *() const
{
    return **my_iter;
}

inline family_base_const_iterator::reference
family_base_const_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline family_base_const_iterator&
family_base_const_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline family_base_const_iterator&
family_base_const_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline family_base_const_iterator&
family_base_const_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline family_base_const_iterator&
family_base_const_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline family_base_const_iterator
family_base_const_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline family_base_const_iterator
family_base_const_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: mate_base_iterator
//============================================================================
//
inline
mate_base_iterator::mate_base_iterator() 
  : my_iter()
{}

inline
mate_base_iterator::mate_base_iterator(const cursor_type& iter)
  : my_iter(iter)
{}

//----------
//
inline bool
mate_base_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
mate_base_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

//----------
//
inline mate_base_iterator::pointer
mate_base_iterator::operator ->() const
{
    return &my_iter->second;
}

inline mate_base_iterator::reference
mate_base_iterator::operator *() const
{
    return my_iter->second;
}

//----------
//
inline mate_base_iterator&
mate_base_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline mate_base_iterator
mate_base_iterator::operator ++(int)
{
    mate_base_iterator  tmp(*this);
    ++my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: mate_base_const_iterator
//============================================================================
//
inline
mate_base_const_iterator::mate_base_const_iterator() 
  : my_iter()
{}

inline
mate_base_const_iterator::mate_base_const_iterator
(const mate_base_iterator& iter)
  : my_iter(iter.my_iter)
{}

inline
mate_base_const_iterator::mate_base_const_iterator(const cursor_type& iter)
  : my_iter(iter)
{}

inline
mate_base_const_iterator::mate_base_const_iterator
(const const_cursor_type& iter)
  : my_iter(iter)
{}

//----------
//
inline bool
mate_base_const_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
mate_base_const_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

//----------
//
inline mate_base_const_iterator::pointer
mate_base_const_iterator::operator ->() const
{
    return &my_iter->second;
}

inline mate_base_const_iterator::reference
mate_base_const_iterator::operator *() const
{
    return my_iter->second;
}

//----------
//
inline mate_base_const_iterator&
mate_base_const_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline mate_base_const_iterator
mate_base_const_iterator::operator ++(int)
{
    mate_base_const_iterator    tmp(*this);
    ++my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: member_base_iterator
//============================================================================
//
inline
member_base_iterator::member_base_iterator() 
  : my_iter()
{}

inline
member_base_iterator::member_base_iterator(const cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
member_base_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
member_base_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
member_base_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline member_base_iterator::pointer
member_base_iterator::operator ->() const
{
    return *my_iter;
}

inline member_base_iterator::reference
member_base_iterator::operator *() const
{
    return **my_iter;
}

inline member_base_iterator::reference
member_base_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline member_base_iterator&
member_base_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline member_base_iterator&
member_base_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline member_base_iterator&
member_base_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline member_base_iterator&
member_base_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline member_base_iterator
member_base_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline member_base_iterator
member_base_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: member_base_const_iterator
//============================================================================
//
inline
member_base_const_iterator::member_base_const_iterator() 
  : my_iter()
{}

inline
member_base_const_iterator::member_base_const_iterator
(const member_base_iterator& iter) 
  : my_iter(iter.my_iter)
{}

inline
member_base_const_iterator::member_base_const_iterator
(const cursor_type& iter) 
  : my_iter(iter)
{}

inline
member_base_const_iterator::member_base_const_iterator
(const const_cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
member_base_const_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
member_base_const_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
member_base_const_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline member_base_const_iterator::pointer
member_base_const_iterator::operator ->() const
{
    return *my_iter;
}

inline member_base_const_iterator::reference
member_base_const_iterator::operator *() const
{
    return **my_iter;
}

inline member_base_const_iterator::reference
member_base_const_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline member_base_const_iterator&
member_base_const_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline member_base_const_iterator&
member_base_const_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline member_base_const_iterator&
member_base_const_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline member_base_const_iterator&
member_base_const_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline member_base_const_iterator
member_base_const_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline member_base_const_iterator
member_base_const_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: offspring_base_iterator
//============================================================================
//
inline
offspring_base_iterator::offspring_base_iterator() 
  : my_data(0)
{}

inline
offspring_base_iterator::offspring_base_iterator(member_id data) 
  : my_data(data)
{}

//----------
//
inline bool
offspring_base_iterator::operator ==(const this_type& c) const
{
    return my_data == c.my_data;
}

inline bool
offspring_base_iterator::operator !=(const this_type& c) const
{
    return my_data != c.my_data;
}

//----------
//
inline offspring_base_iterator::pointer
offspring_base_iterator::operator ->() const
{
    return my_data;
}

inline offspring_base_iterator::reference
offspring_base_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline offspring_base_iterator&
offspring_base_iterator::operator ++()
{
    advance();
    return *this;
}

inline offspring_base_iterator
offspring_base_iterator::operator ++(int)
{
    offspring_base_iterator tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: offspring_base_const_iterator
//============================================================================
//
inline
offspring_base_const_iterator::offspring_base_const_iterator() 
  : my_data(0)
{}

inline
offspring_base_const_iterator::offspring_base_const_iterator
(const offspring_base_iterator& iter) 
  : my_data(iter.my_data)
{}

inline
offspring_base_const_iterator::offspring_base_const_iterator(member_id data) 
  : my_data(data)
{}

//----------
//
inline bool
offspring_base_const_iterator::operator ==(const this_type& c) const
{
    return my_data == c.my_data;
}

inline bool
offspring_base_const_iterator::operator !=(const this_type& c) const
{
    return my_data != c.my_data;
}

//----------
//
inline offspring_base_const_iterator::pointer
offspring_base_const_iterator::operator ->() const
{
    return my_data;
}

inline offspring_base_const_iterator::reference
offspring_base_const_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline offspring_base_const_iterator&
offspring_base_const_iterator::operator ++()
{
    advance();
    return *this;
}

inline offspring_base_const_iterator
offspring_base_const_iterator::operator ++(int)
{
    offspring_base_const_iterator   tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: parent_base_iterator
//============================================================================
//
inline
parent_base_iterator::parent_base_iterator() 
  : my_family(0), my_data(0)
{}

inline
parent_base_iterator::parent_base_iterator(family_id data) 
  : my_family(data), my_data(0)
{
    init(); 
}

//----------
//
inline bool
parent_base_iterator::operator ==(const this_type& c) const
{
    return my_data == c.my_data  &&  my_family == c.my_family;
}

inline bool
parent_base_iterator::operator !=(const this_type& c) const
{
    return my_data != c.my_data  ||  my_family != c.my_family;
}

//----------
//
inline parent_base_iterator::pointer
parent_base_iterator::operator ->() const
{
    return my_data;
}

inline parent_base_iterator::reference
parent_base_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline parent_base_iterator&
parent_base_iterator::operator ++()
{
    advance();
    return *this;
}

inline parent_base_iterator
parent_base_iterator::operator ++(int)
{
    parent_base_iterator    tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: parent_base_const_iterator
//============================================================================
//
inline
parent_base_const_iterator::parent_base_const_iterator() 
  : my_family(0), my_data(0)
{}

inline
parent_base_const_iterator::parent_base_const_iterator
(const parent_base_iterator& iter) 
  : my_family(iter.my_family), my_data(iter.my_data)
{}

inline
parent_base_const_iterator::parent_base_const_iterator(family_id data) 
  : my_family(data), my_data(0)
{
    init(); 
}

//----------
//
inline bool
parent_base_const_iterator::operator ==(const this_type& c) const
{
    return my_data == c.my_data  &&  my_family == c.my_family;
}

inline bool
parent_base_const_iterator::operator !=(const this_type& c) const
{
    return my_data != c.my_data  ||  my_family != c.my_family;
}

//----------
//
inline parent_base_const_iterator::pointer
parent_base_const_iterator::operator ->() const
{
    return my_data;
}

inline parent_base_const_iterator::reference
parent_base_const_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline parent_base_const_iterator&
parent_base_const_iterator::operator ++()
{
    advance();
    return *this;
}

inline parent_base_const_iterator
parent_base_const_iterator::operator ++(int)
{
    parent_base_const_iterator  tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: pedigree_base_iterator
//============================================================================
//
inline
pedigree_base_iterator::pedigree_base_iterator() 
  : my_iter()
{}

inline
pedigree_base_iterator::pedigree_base_iterator(const cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
pedigree_base_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
pedigree_base_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
pedigree_base_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline pedigree_base_iterator::pointer
pedigree_base_iterator::operator ->() const
{
    return *my_iter;
}

inline pedigree_base_iterator::reference
pedigree_base_iterator::operator *() const
{
    return **my_iter;
}

inline pedigree_base_iterator::reference
pedigree_base_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline pedigree_base_iterator&
pedigree_base_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline pedigree_base_iterator&
pedigree_base_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline pedigree_base_iterator&
pedigree_base_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline pedigree_base_iterator&
pedigree_base_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline pedigree_base_iterator
pedigree_base_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline pedigree_base_iterator
pedigree_base_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: pedigree_base_const_iterator
//============================================================================
//
inline
pedigree_base_const_iterator::pedigree_base_const_iterator() 
  : my_iter()
{}

inline
pedigree_base_const_iterator::pedigree_base_const_iterator
(const pedigree_base_iterator& iter) 
  : my_iter(iter.my_iter)
{}

inline
pedigree_base_const_iterator::pedigree_base_const_iterator
(const cursor_type& iter) 
  : my_iter(iter)
{}

inline
pedigree_base_const_iterator::pedigree_base_const_iterator
(const const_cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
pedigree_base_const_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
pedigree_base_const_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
pedigree_base_const_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline pedigree_base_const_iterator::pointer
pedigree_base_const_iterator::operator ->() const
{
    return *my_iter;
}

inline pedigree_base_const_iterator::reference
pedigree_base_const_iterator::operator *() const
{
    return **my_iter;
}

inline pedigree_base_const_iterator::reference
pedigree_base_const_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline pedigree_base_const_iterator&
pedigree_base_const_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline pedigree_base_const_iterator&
pedigree_base_const_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline pedigree_base_const_iterator&
pedigree_base_const_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline pedigree_base_const_iterator&
pedigree_base_const_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline pedigree_base_const_iterator
pedigree_base_const_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline pedigree_base_const_iterator
pedigree_base_const_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: progeny_base_iterator
//============================================================================
//
inline
progeny_base_iterator::progeny_base_iterator() 
  : my_data(0), my_iter(), my_limit()
{}

inline 
progeny_base_iterator::progeny_base_iterator
(const cursor_type& iter, const cursor_type& limit)
  : my_data(0), my_iter(iter), my_limit(limit)
{
    init();
}

//----------
//
inline bool
progeny_base_iterator::operator ==(const progeny_base_iterator& c) const
{
    return  my_data == c.my_data  &&  my_iter == c.my_iter  &&
            my_limit == c.my_limit;
}

inline bool
progeny_base_iterator::operator !=(const progeny_base_iterator& c) const
{
    return  my_data != c.my_data  ||  my_iter != c.my_iter  ||
            my_limit != c.my_limit;
}

//----------
//
inline progeny_base_iterator::pointer
progeny_base_iterator::operator ->() const
{
    return my_data;
}

inline progeny_base_iterator::reference
progeny_base_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline progeny_base_iterator&
progeny_base_iterator::operator ++()
{
    advance();
    return *this;
}

inline progeny_base_iterator
progeny_base_iterator::operator ++(int)
{
    progeny_base_iterator   tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: progeny_base_const_iterator
//============================================================================
//
inline
progeny_base_const_iterator::progeny_base_const_iterator() 
  : my_data(0), my_iter(), my_limit()
{}

inline 
progeny_base_const_iterator::progeny_base_const_iterator
(const progeny_base_iterator& iter)
  : my_data(0), my_iter(iter.my_iter), my_limit(iter.my_limit)
{}

inline 
progeny_base_const_iterator::progeny_base_const_iterator
(const cursor_type& iter, const cursor_type& limit)
  : my_data(0), my_iter(iter), my_limit(limit)
{
    init();
}

inline 
progeny_base_const_iterator::progeny_base_const_iterator
(const const_cursor_type& iter, const const_cursor_type& limit)
  : my_data(0), my_iter(iter), my_limit(limit)
{
    init();
}

//----------
//
inline bool
progeny_base_const_iterator::operator ==(const this_type& c) const
{
    return  my_data == c.my_data  &&  my_iter == c.my_iter  &&
            my_limit == c.my_limit;
}

inline bool
progeny_base_const_iterator::operator !=(const this_type& c) const
{
    return  my_data != c.my_data  ||  my_iter != c.my_iter  ||
            my_limit != c.my_limit;
}

//----------
//
inline progeny_base_const_iterator::pointer
progeny_base_const_iterator::operator ->() const
{
    return my_data;
}

inline progeny_base_const_iterator::reference
progeny_base_const_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline progeny_base_const_iterator&
progeny_base_const_iterator::operator ++()
{
    advance();
    return *this;
}

inline progeny_base_const_iterator
progeny_base_const_iterator::operator ++(int)
{
    progeny_base_const_iterator tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: sibling_base_iterator
//============================================================================
//
inline
sibling_base_iterator::sibling_base_iterator() 
  : my_self(0), my_data(0)
{}

inline
sibling_base_iterator::sibling_base_iterator(member_id self) 
  : my_self(self), my_data(0)
{}

inline
sibling_base_iterator::sibling_base_iterator(member_id self, member_id data) 
  : my_self(self), my_data(data)
{
    if (my_self  &&  my_data == my_self)
    {
        advance();
    }
}

//----------
//
inline bool
sibling_base_iterator::operator ==(const this_type& c) const
{
    return my_self == c.my_self  &&  my_data == c.my_data;
}

inline bool
sibling_base_iterator::operator !=(const this_type& c) const
{
    return my_self != c.my_self  ||  my_data != c.my_data;
}

//----------
//
inline sibling_base_iterator::pointer
sibling_base_iterator::operator ->() const
{
    return my_data;
}

inline sibling_base_iterator::reference
sibling_base_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline sibling_base_iterator&
sibling_base_iterator::operator ++()
{
    advance();
    return *this;
}

inline sibling_base_iterator
sibling_base_iterator::operator ++(int)
{
    sibling_base_iterator   tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: sibling_base_const_iterator
//============================================================================
//
inline
sibling_base_const_iterator::sibling_base_const_iterator() 
  : my_self(0), my_data(0)
{}

inline
sibling_base_const_iterator::sibling_base_const_iterator
(const sibling_base_iterator& iter) 
  : my_self(iter.my_self), my_data(iter.my_data)
{}

inline
sibling_base_const_iterator::sibling_base_const_iterator(member_id self) 
  : my_self(self), my_data(0)
{}

inline
sibling_base_const_iterator::sibling_base_const_iterator
(member_id self, member_id data) 
  : my_self(self), my_data(data)
{
    if (my_self  &&  my_data == my_self)
    {
        advance();
    }
}

//----------
//
inline bool
sibling_base_const_iterator::operator ==(const this_type& c) const
{
    return my_self == c.my_self  &&  my_data == c.my_data;
}

inline bool
sibling_base_const_iterator::operator !=(const this_type& c) const
{
    return my_self != c.my_self  ||  my_data != c.my_data;
}

//----------
//
inline sibling_base_const_iterator::pointer
sibling_base_const_iterator::operator ->() const
{
    return my_data;
}

inline sibling_base_const_iterator::reference
sibling_base_const_iterator::operator *() const
{
    return *my_data;
}

//----------
//
inline sibling_base_const_iterator&
sibling_base_const_iterator::operator ++()
{
    advance();
    return *this;
}

inline sibling_base_const_iterator
sibling_base_const_iterator::operator ++(int)
{
    sibling_base_const_iterator tmp(*this);
    advance();
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: subpedigree_base_iterator
//============================================================================
//
inline
subpedigree_base_iterator::subpedigree_base_iterator() 
  : my_iter()
{}

inline
subpedigree_base_iterator::subpedigree_base_iterator(const cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
subpedigree_base_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
subpedigree_base_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
subpedigree_base_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline subpedigree_base_iterator::pointer
subpedigree_base_iterator::operator ->() const
{
    return *my_iter;
}

inline subpedigree_base_iterator::reference
subpedigree_base_iterator::operator *() const
{
    return **my_iter;
}

inline subpedigree_base_iterator::reference
subpedigree_base_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline subpedigree_base_iterator&
subpedigree_base_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline subpedigree_base_iterator&
subpedigree_base_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline subpedigree_base_iterator&
subpedigree_base_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline subpedigree_base_iterator&
subpedigree_base_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline subpedigree_base_iterator
subpedigree_base_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline subpedigree_base_iterator
subpedigree_base_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}


//============================================================================
//  IMPLEMENTATION: subpedigree_base_const_iterator
//============================================================================
//
inline
subpedigree_base_const_iterator::subpedigree_base_const_iterator() 
  : my_iter()
{}

inline
subpedigree_base_const_iterator::subpedigree_base_const_iterator
(const subpedigree_base_iterator& iter) 
  : my_iter(iter.my_iter)
{}

inline
subpedigree_base_const_iterator::subpedigree_base_const_iterator
(const cursor_type& iter) 
  : my_iter(iter)
{}

inline
subpedigree_base_const_iterator::subpedigree_base_const_iterator
(const const_cursor_type& iter) 
  : my_iter(iter)
{}

//----------
//
inline bool
subpedigree_base_const_iterator::operator ==(const this_type& c) const
{
    return my_iter == c.my_iter;
}

inline bool
subpedigree_base_const_iterator::operator !=(const this_type& c) const
{
    return my_iter != c.my_iter;
}

inline bool
subpedigree_base_const_iterator::operator <(const this_type& c) const
{
    return my_iter < c.my_iter;
}

//----------
//
inline subpedigree_base_const_iterator::pointer
subpedigree_base_const_iterator::operator ->() const
{
    return *my_iter;
}

inline subpedigree_base_const_iterator::reference
subpedigree_base_const_iterator::operator *() const
{
    return **my_iter;
}

inline subpedigree_base_const_iterator::reference
subpedigree_base_const_iterator::operator [](difference_type n) const
{
    return *(my_iter[n]);
}

//----------
//
inline subpedigree_base_const_iterator&
subpedigree_base_const_iterator::operator +=(difference_type n)
{
    my_iter += n;
    return *this;
}

inline subpedigree_base_const_iterator&
subpedigree_base_const_iterator::operator -=(difference_type n)
{
    my_iter -= n;
    return *this;
}

inline subpedigree_base_const_iterator&
subpedigree_base_const_iterator::operator ++()
{
    ++my_iter;
    return *this;
}

inline subpedigree_base_const_iterator&
subpedigree_base_const_iterator::operator --()
{
    --my_iter;
    return *this;
}

inline subpedigree_base_const_iterator
subpedigree_base_const_iterator::operator ++(int)
{
    this_type   tmp(*this);
    ++my_iter;
    return tmp;
}

inline subpedigree_base_const_iterator
subpedigree_base_const_iterator::operator --(int)
{
    this_type   tmp(*this);
    --my_iter;
    return tmp;
}

} // End namespace MPED
} // End namespace SAGE
