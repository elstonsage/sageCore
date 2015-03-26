#ifndef __BIT_FIELD_H
#define __BIT_FIELD_H

//
//  Bit field 1.0 - similar to the bit vector, without insertion or
//                  resizing, but adding set operations, binary operations (&,
//                  &=, |, etc) and efficient comparisions.
//
//  Author:  Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History: 0.1 gcw Initial Implementation               Aug    1997
//           1.0 gcw Revised and Finished                 Sep  4 1997
//           1.1 gcw Got rid of vector<bool> dependency   Nov 17 1997
//           1.2 gcw Added many features that should have Mar 23 1998
//                   been there in the first place
//           1.3 kbj Optimized several of the operators   Mar 30 1998
//
//  Note: Do not assume that vector operations (empty, insert, etc) work.
//        Notably empty(), clear,  and the insertion/deletion operations work
//        differently or are not available.
//
//  
//  Set Operations: 
//    Set operations Union, Intersection and Difference modify the bit_field
//    evaluated.  Non-modifying versions are included as non-member
//    functions. The modifying and non-modifying versions of most of the
//    optimized bit functions are as follows:
//
//    Modifying           Non-modifying
//    ---------           -------------
//    a &= b              a & b
//    a |= b              a | b
//    a ^= b              a ^ b
//    a.lshift(n)         a << n
//    a.rshift(n)         a >> n
//    a.Intersection(b)   Intersection(a, b)  (equiv. to  a &  b)
//    a.Union(b)          Union(a, b)         (equiv. to  a |  b)
//    a.Difference(b)     Difference(a, b)    (equiv. to  a & ~b)
//    a.flip(b)           ~a
//    a.And(b)            And(a, b)
//    a.Or(b)             Or(a, b)
//    a.Xor(b)            Xor(a, b)
//    a += b              a + b               (addition as if field was number)
//    a -= b              a - b               (subtraction "")
//
//                        a.subset(b)
//                        a.disjoint(b)
//                        a.complement(b)
//                        a.empty()
//
//    The functions And, Or and Xor (member and non-member versions) provide
//    specialized versions of operators (&, |, ^, etc).  In the operator
//    based versions, the size of the resulting bit_field is equal to the
//    size of the smaller of the two bit_field inputs.  In the specialized
//    versions, the size of the result is equal to that of the first
//    element.  The second element, if smaller, is assumed to contain zeros in the
//    digits not assigned.
//
//  Copyright (c) 1998  R.C. Elston

#include <string>
#include <vector>
#include <string.h>
#include <utility>
#include <algorithm>

using namespace std;

// __U32 is an unsigned 32 bit integer.
#if __KCC || __DECCXX
#define __U32 unsigned int
#else
#define __U32 unsigned int
#endif

#define __WORD_BIT (int(CHAR_BIT*sizeof(__U32)))

class bit_field_iterator;
class const_bit_field_iterator;

class bit_field_reference
{
public:

  friend class bit_field_iterator;
  friend class const_bit_field_iterator;

  typedef bit_field_reference reference;

// Constructor

  bit_field_reference() : p(0), mask(0) {}

// Conversion, etc.

  operator bool() const { return !!(*p & mask); }

  reference& operator=(bool x)
  { if (x) *p |= mask;
    else   *p &= ~mask;
    return *this;
  }
  reference& operator=(const reference& x) { return *this = bool(x); }

// Comparison Operators

  bool operator==(const reference& x) const { return bool(*this)== bool(x); }
  bool operator< (const reference&  ) const { return false; }

  void flip() { *p ^= mask; }

protected:

  unsigned int* p;
  unsigned int mask;

// Private Constructor

  bit_field_reference(unsigned int* x, unsigned int y) : p(x), mask(y) {}
};


class bit_field_iterator
{
public:

  friend class bit_field;
  friend class const_bit_field_iterator;

  typedef bit_field_iterator       iterator;
  typedef const_bit_field_iterator const_iterator;
  typedef bit_field_reference      reference;
  typedef bit_field_reference*     pointer;
  typedef bool                     const_reference;
  typedef bool                     value_type;
  typedef ptrdiff_t                difference_type;
  typedef size_t                   size_type;
  typedef bidirectional_iterator_tag iterator_category;

// Constructors

  bit_field_iterator() : p(0), offset(0) {}
  bit_field_iterator(const bit_field_iterator& b) : p(b.p), offset(b.offset) {}
  bit_field_iterator(unsigned int* x, unsigned int y) : p(x), offset(y) {}

// Dereference

  inline reference operator*() const;

  inline reference operator[](difference_type i) const;

// Increment/Decrement

  inline iterator& operator++();
  inline iterator  operator++(int);
  inline iterator& operator--();
  inline iterator  operator--(int);

  inline iterator& operator+=(difference_type i);
  inline iterator& operator-=(difference_type i);

  inline iterator operator+(difference_type i) const; 
  inline iterator operator-(difference_type i) const;

  inline difference_type operator-(iterator x) const;

// Comparison

  inline bool operator==(const iterator& x) const;
  inline bool operator!=(const iterator& x) const;
  
  inline bool operator<(iterator x) const;

protected:

  inline void bump_up();
  inline void bump_down();

  unsigned int* p;
  unsigned int  offset;

};

class const_bit_field_iterator
{
public:

  friend class bit_field;

  typedef bit_field_iterator       iterator;
  typedef const_bit_field_iterator const_iterator;
  typedef bit_field_reference      reference;
  typedef bit_field_reference*     pointer;
  typedef bool                     value_type;
  typedef bool                     const_reference;
  typedef ptrdiff_t                difference_type;
  typedef size_t                   size_type;
  typedef bidirectional_iterator_tag iterator_category;
// Constructors

  const_bit_field_iterator() : p(0), offset(0) {}
  const_bit_field_iterator(unsigned int* x, unsigned int y) : p(x), offset(y) {}
  const_bit_field_iterator(const iterator& x) : p(x.p), offset(x.offset) {}

// Dereferencing

  inline const_reference operator*() const;

  const_reference operator[](difference_type i) const;

// Increment/Decrement

  inline const_iterator& operator++();
  inline const_iterator  operator++(int);
  inline const_iterator& operator--();
  inline const_iterator  operator--(int); 

  inline const_iterator& operator+=(difference_type i);
  inline const_iterator& operator-=(difference_type i);

  inline const_iterator operator+(difference_type i) const;
  inline const_iterator operator-(difference_type i) const;

  inline difference_type operator-(const_iterator x) const;

// Comparisons

  inline bool operator==(const const_iterator& x) const;
  inline bool operator!=(const const_iterator& x) const;

  inline bool operator<(const_iterator x) const;

protected:

  inline void bump_up();
  inline void bump_down();

  unsigned int* p;
  unsigned int  offset;

};

class bit_field
{
public:

  typedef __U32                    mem_type;
  typedef bit_field_reference      reference;
  typedef bool                     const_reference;
  typedef bool                     value_type;
  typedef ptrdiff_t                difference_type;
  typedef size_t                   size_type;  
  typedef bit_field_iterator       iterator;
  typedef const_bit_field_iterator const_iterator;

// Constructors/Destructor
  bit_field() { /*cout << "allocating 1: " << this << endl;*/ initialize(0); }
  explicit bit_field(size_type n, bool value = bool() )
  { /*cout << "allocating 2: " << this << endl;*/ initialize(n); ::memset(mbegin(),(value? ~0:0), bytesize()); }
  explicit bit_field(size_type n, __U32 u)
  {
    /*cout << "allocating 3: " << this << endl;*/
    initialize(n);
    ::memset(mbegin(),0, bytesize());
    if(n && u) *mbegin()       = u;
  }
  bit_field(const std::vector<bool>& b)
  { /*cout << "allocating 4: " << this << endl;*/ initialize(b.size()); std::copy(b.begin(), b.end(), begin()); }
  bit_field(const bit_field& b)
  { /*cout << "allocating 5: " << this << endl;*/ initialize(b.size()); ::memcpy(mbegin(),b.mbegin(), bytesize()); }
  bit_field(const_iterator s, const_iterator e)
  {
    /*cout << "allocating 6: " << this << endl;*/ 
    initialize(e-s);
    if(s.offset == 0)
      ::memcpy(mbegin(),s.p, bytesize()); 
    else
    {
      iterator i = begin();
      for( ; s != e; ++s, ++i)
        *i = *s;
    }
  }
  bit_field(const void* s, const void* e)
  { 
    /*cout << "allocating 7: " << this << endl;*/ 
    size_t n = (char*)e-(char*)s;
    initialize(n*CHAR_BIT); 
    memcpy(mbegin(),s,n); 
  }
  bit_field(const bool* s, const bool* e)
  { /*cout << "allocating 8: " << this << endl;*/ initialize(e - s); std::copy(s, e, begin()); }

  ~bit_field() { /*cout << "deallocating: " << this << ' ' << start.p << endl;*/ deallocate(); }

  inline bit_field& operator= (const bit_field& rhs);

// Iterators and references

  inline reference       operator[] (size_type n);
  inline const_reference operator[] (size_type n) const;

  inline iterator       begin();
  inline const_iterator begin() const;
  inline iterator         end();
  inline const_iterator   end() const;
  
  inline reference       front();
  inline const_reference front() const;
  inline reference        back();
  inline const_reference  back() const;
  
// Miscellaneous operations

  size_type     size() const;
  size_type max_size() const;
  size_type capacity() const;

  // NOTE: empty returns if any bits are set, so is with set operations!!

  void swap(bit_field& rhs);
  
// Conversion operations

  __U32  to_u16(size_type t = 0) const;  // masked to 16 bits
  __U32  to_u32(size_type t = 0) const;
  string to_string() const;

  __U32  chunk(size_type t = 0) const;

  void or_chunk(__U32 u, unsigned int offset);

  size_type chunk_size() const;

  void set_chunk(__U32, size_type t = 0);
  void set_u32(__U32, size_type t = 0);
  void set_u16(__U32, size_type t = 0);

// Operators

  // Equality operators

  inline bool operator== (const bit_field& rhs) const;
  inline bool operator!= (const bit_field& rhs) const;

  // Increment/Decrement Operators

  inline bit_field& operator++();
  inline bit_field  operator++(int);
  inline bit_field& operator--();
  inline bit_field  operator--(int);
  
  inline bit_field  operator+  (const bit_field& rhs) const;
  inline bit_field& operator+= (const bit_field& rhs);
  inline bit_field  operator-  (const bit_field& rhs) const;
  inline bit_field& operator-= (const bit_field& rhs);

  // Bitwise operators

  inline bit_field  operator&  (const bit_field& rhs) const;
  inline bit_field& operator&= (const bit_field& rhs);
  inline bit_field  operator|  (const bit_field& rhs) const;
  inline bit_field& operator|= (const bit_field& rhs);
  inline bit_field  operator^  (const bit_field& rhs) const;
  inline bit_field& operator^= (const bit_field& rhs);
  
  inline bit_field operator<< (size_type t) const;
  inline bit_field operator>> (size_type t) const;
  
  inline bit_field& operator<<= (size_type t);
  inline bit_field& operator>>= (size_type t);

  inline bit_field& lshift(size_type t);
  inline bit_field& rshift(size_type t);
 
  inline bit_field operator~ () const;

  inline bit_field& flip();
  inline bit_field& clear();

// Set operations

  inline bit_field& And          (const bit_field& rhs);
  inline bit_field& Or           (const bit_field& rhs);
  inline bit_field& Xor          (const bit_field& rhs);
  inline bit_field& Union        (const bit_field& rhs);
  inline bit_field& Intersection (const bit_field& rhs);
  inline bit_field& Difference   (const bit_field& rhs);

  inline bool subset    (const bit_field& rhs) const;
  inline bool disjoint  (const bit_field& rhs) const;
  inline bool complement(const bit_field& rhs) const;

  inline bool empty() const;
  inline bool any  () const;
  
// Memory allocation, etc

  inline void reserve(size_type t);

  inline void resize(size_type t, bool b = false)
  {
    if(t <= (unsigned) (finish - start)) { finish = start + t; return; }
    reserve(t);
    if(!b)
    { set_zero(); for(mem_type* m = mfinish()+1; m < mend(); ++m) *m = 0; }
    else
    { set_one (); for(mem_type* m = mfinish()+1; m < mend(); ++m) *m = (mem_type) -1; }
    finish = start + t;
  } 
  
protected:

  // Note:  Data allocation schemes vary wildly for different compilers, so
  // no real 'allocator' is currently used.  This may change when compiler
  // support of the standard STL allocators improves.  See Implementation
  // for more details
  inline unsigned int* bit_alloc(size_type n);

  inline void deallocate();

  inline void initialize(size_type n);

// Internal plus, minus, and, or, and xor

  inline void int_plus  (mem_type* b1, mem_type* e1, mem_type* b2);
  inline void int_minus (mem_type* b1, mem_type* e1, mem_type* b2);
  inline void int_and   (mem_type* b1, mem_type* e1, mem_type* b2);
  inline void int_or    (mem_type* b1, mem_type* e1, mem_type* b2);
  inline void int_xor   (mem_type* b1, mem_type* e1, mem_type* b2);

// Memory Member functions to make things easier.

  inline mem_type*  mbegin() const;
  inline mem_type*    mend() const;
  inline mem_type* mfinish() const;
  
  inline size_type bytesize() const;

  inline mem_type get_mask() const;

// Set all past end to zero or one

  void set_zero() const;
  void set_one()  const;

// Data Elements

  iterator start;
  iterator finish;

  unsigned int* end_of_storage;
};

// ================
// Inline Functions
// ================

inline bit_field Union(const bit_field& lhs, const bit_field& rhs)
{
  return lhs | rhs;
}

inline bit_field Intersection(const bit_field& lhs, const bit_field& rhs)
{
  return lhs & rhs;
}

inline bit_field Difference(const bit_field& lhs, const bit_field& rhs)
{
  return lhs & ~rhs;
}

inline bit_field And(const bit_field& lhs, const bit_field& rhs)
{ bit_field b = lhs; return b.And(rhs); }

inline bit_field Or(const bit_field& lhs, const bit_field& rhs)
{ bit_field b = lhs; return b.Or(rhs); }

inline bit_field Xor(const bit_field& lhs, const bit_field& rhs)
{ bit_field b = lhs; return b.Xor(rhs); }

// =======================
// Inline Member Functions
// =======================

// ==================
// bit_field_iterator
// ==================

inline bit_field_iterator::reference bit_field_iterator::operator*() const
{ return reference(p, 1U << offset); }

inline bit_field_iterator::reference
    bit_field_iterator::operator[](difference_type i) const
{ return *(*this + i); }

inline bit_field_iterator& bit_field_iterator::operator++()
{ bump_up(); return *this; }

inline bit_field_iterator  bit_field_iterator::operator++(int)
{ iterator tmp = *this; bump_up(); return tmp; }

inline bit_field_iterator& bit_field_iterator::operator--()
{ bump_down(); return *this; }

inline bit_field_iterator  bit_field_iterator::operator--(int)
{ iterator tmp = *this; bump_down(); return tmp; }

inline bit_field_iterator& bit_field_iterator::operator+=(difference_type i)
{
  difference_type n = i + offset;
  p += n / __WORD_BIT;
  n = n % __WORD_BIT;
  if (n < 0) { offset = n + __WORD_BIT; --p; }
  else       { offset = n;                   }
  return *this;
}

inline bit_field_iterator& bit_field_iterator::operator-=(difference_type i)
{ *this += -i; return *this; }

inline bit_field_iterator bit_field_iterator::operator+(difference_type i) const 
{ iterator tmp = *this; return tmp += i; }

inline bit_field_iterator bit_field_iterator::operator-(difference_type i) const
{ iterator tmp = *this; return tmp -= i; }

inline bit_field_iterator::difference_type
    bit_field_iterator::operator-(iterator x) const
{ return __WORD_BIT * (p - x.p) + offset - x.offset; }

inline bool bit_field_iterator::operator==(const iterator& x) const 
{ return p == x.p && offset == x.offset; }

inline bool bit_field_iterator::operator!=(const iterator& x) const
{ return p != x.p || offset != x.offset; }

inline bool bit_field_iterator::operator<(iterator x) const 
{ return p < x.p || (p == x.p && offset < x.offset); }

inline void bit_field_iterator::bump_up()
{ if (offset++ == __WORD_BIT - 1) { offset = 0; ++p; } }

inline void bit_field_iterator::bump_down()
{ if (offset-- == 0) { offset = __WORD_BIT - 1; --p; } }

// ========================
// const_bit_field_iterator
// ========================

inline const_bit_field_iterator::const_reference
     const_bit_field_iterator::operator*() const
{ return reference(p, 1U << offset); }

inline const_bit_field_iterator::const_reference
    const_bit_field_iterator::operator[](difference_type i) const 
{ return *(*this + i); }

inline const_bit_field_iterator& const_bit_field_iterator::operator++()
{ bump_up(); return *this; }

inline const_bit_field_iterator const_bit_field_iterator::operator++(int)
{ const_iterator tmp = *this; bump_up(); return tmp; }

inline const_bit_field_iterator& const_bit_field_iterator::operator--()
{ bump_down(); return *this; }

inline const_bit_field_iterator const_bit_field_iterator::operator--(int) 
{ const_iterator tmp = *this; bump_down(); return tmp; }

inline const_bit_field_iterator&
    const_bit_field_iterator::operator+=(difference_type i)
{
  difference_type n = i + offset;
  p += n / __WORD_BIT;
  n = n % __WORD_BIT;
  if (n < 0) { offset = n + __WORD_BIT; --p; }
  else       { offset = n;                   }
  return *this;
}

inline const_bit_field_iterator&
    const_bit_field_iterator::operator-=(difference_type i)
{ *this += -i; return *this; }

inline const_bit_field_iterator const_bit_field_iterator::operator+(difference_type i) const
{ const_iterator tmp = *this; return tmp += i; }

inline const_bit_field_iterator const_bit_field_iterator::operator-(difference_type i) const
{ const_iterator tmp = *this; return tmp -= i; }

inline const_bit_field_iterator::difference_type
    const_bit_field_iterator::operator-(const_iterator x) const 
{ return __WORD_BIT * (p - x.p) + offset - x.offset; }

inline bool const_bit_field_iterator::operator==(const const_iterator& x) const 
{ return p == x.p && offset == x.offset; }

inline bool const_bit_field_iterator::operator!=(const const_iterator& x) const
{ return p != x.p || offset != x.offset; }

inline bool const_bit_field_iterator::operator<(const_iterator x) const
{ return p < x.p || (p == x.p && offset < x.offset); }

inline void const_bit_field_iterator::bump_up()   { if (offset++ == __WORD_BIT - 1) { offset = 0; ++p; } }
inline void const_bit_field_iterator::bump_down() { if (offset-- == 0) { offset = __WORD_BIT - 1; --p; } }

// ==========================
// bit_field Inline Functions
// ==========================

inline bit_field& bit_field::operator= (const bit_field& rhs)
{
  if(this == &rhs) return *this;
  if(rhs.size() > capacity())
  {
    deallocate();
    initialize(rhs.size());
  }
  else
    finish = begin() + rhs.size();
  ::memcpy(mbegin(),rhs.mbegin(), bytesize());
  return *this;
}

inline bit_field::reference bit_field::operator[] (size_type n)
{ return *(begin() + n); }

inline bit_field::const_reference bit_field::operator[] (size_type n) const
{ return *(begin() + n); }

inline bit_field::iterator       bit_field::begin()       { return start;  }
inline bit_field::const_iterator bit_field::begin() const { return start;  }
inline bit_field::iterator       bit_field::end  ()       { return finish; }
inline bit_field::const_iterator bit_field::end  () const { return finish; }

inline bit_field::reference       bit_field::front()       { return *begin();     }
inline bit_field::const_reference bit_field::front() const { return *begin();     }
inline bit_field::reference       bit_field::back ()       { return *(end() - 1); }
inline bit_field::const_reference bit_field::back () const { return *(end() - 1); }

inline bit_field::size_type bit_field::size() const
{ return size_type(end()-begin()); }

inline bit_field::size_type bit_field::max_size() const
{ return size_type(-1); }

inline bit_field::size_type bit_field::capacity() const 
{ return size_type(const_iterator(end_of_storage, 0) - begin()); }

inline void bit_field::swap(bit_field& rhs)
{ 
  std::swap(start, rhs.start);
  std::swap(finish, rhs.finish);
  std::swap(end_of_storage, rhs.end_of_storage);
}

inline __U32 bit_field::to_u16(size_type t) const
{
  __U32 u;
  const_iterator i = begin() + t;
  if(i.p >= mend()) return 0;

  set_zero();

  if(i.offset <= 16)
  {
    u = *i.p >> i.offset;

    return u & ((1 << 16) - 1);
  }
  
  u = (*i.p) >> i.offset;
  
  ++i.p;
  
  if(i.p == mend()) return u;
  
  u |= ((__U32) *i.p) << (16 - i.offset);

  u &= (1 << 16) - 1;
  
  return u;
}

inline __U32 bit_field::to_u32(size_type t) const
{
  __U32 u;
  const_iterator i = begin() + t;
  if(i.p >= mend()) return 0;

  set_zero();

  if(i.offset == 0)
  {
    u = *i.p;
    return u;
  }
  
  u = (*i.p) >> i.offset;
  
  ++i.p;
  
  if(i.p == mend()) return u;
  
  u |= ((__U32) *i.p) << (__WORD_BIT - i.offset);
  
  return u;
}

inline void bit_field::set_u32(__U32 u, size_type t)
{
  iterator i = begin() + t;
  if(i.p >= mend()) return;

  if(i.offset == 0)
  {
    *i.p = u;
    return;
  }

  __U32 umask = 1 << i.offset - 1;
  
  *i.p = (*i.p & umask) | u << i.offset;
  
  ++i.p;
  
  if(i.p == mend()) return;

  umask = ~umask;
  
  *i.p = (*i.p & umask) | u >> (__WORD_BIT - i.offset);
}

inline void bit_field::set_chunk(__U32 u, size_type t)
{
  set_u32(u,t);
}

inline string bit_field::to_string() const
{
  string s(size(), '1');
  for(size_t i = 0; i < size(); ++i)
    if((*this)[i] == 0) s[size() - i - 1] = '0';
  return s;
}

inline __U32 bit_field::chunk(size_type t) const
{
  return to_u32(t);
}

inline bit_field::size_type bit_field::chunk_size() const
{
  return __WORD_BIT;
}

inline bool bit_field::operator== (const bit_field& rhs) const
{
  if(&rhs == this) return true;
  if(size() != rhs.size()) return false;
  set_zero(); rhs.set_zero();
  
  return (::memcmp(mbegin(), rhs.mbegin(), bytesize()) == 0);
}

inline bool bit_field::operator!= (const bit_field& rhs) const
{ return !(*this == rhs); }

inline bit_field& bit_field::operator++()
{
  mem_type* m = mbegin();
  for( ; m != mend() && (*m) + 1 == 0; ++m)
    ++(*m);
  
  if(m != mend()) ++(*m);
    
  return *this;
}

inline bit_field bit_field::operator++(int)
{ bit_field temp = *this; ++*this; return temp; }

inline bit_field& bit_field::operator--()
{
  mem_type* m = mbegin();
  for( ; m != mend() && *m == 0; ++m)
    --(*m);
  
  if(m != mend()) --(*m);
    
  return *this;
}

inline bit_field bit_field::operator--(int)
{ bit_field temp = *this; --*this; return temp; }

inline bit_field bit_field::operator+ (const bit_field& rhs) const
{ bit_field b = *this;  return b += rhs; }

inline bit_field& bit_field::operator+= (const bit_field& rhs)
{ 
  if(size() > rhs.size()) finish = begin() + rhs.size();
  int_plus(mbegin(), mend(), rhs.mbegin());
  return *this;
}

inline bit_field bit_field::operator- (const bit_field& rhs) const
{ bit_field b = *this;  return b -= rhs; }

inline bit_field& bit_field::operator-= (const bit_field& rhs)
{ 
  if(size() > rhs.size()) finish = begin() + rhs.size();
  int_minus(mbegin(), mend(), rhs.mbegin());
  return *this;
}

inline bit_field bit_field::operator& (const bit_field& rhs) const
{ bit_field b = *this;
  return b &= rhs; }

inline bit_field& bit_field::operator&= (const bit_field& rhs)
{
  if(&rhs == this) return *this;
  if(size() > rhs.size()) finish = begin() + rhs.size();
  int_and(mbegin(), mend(), rhs.mbegin());
  return *this;
}

inline bit_field bit_field::operator| (const bit_field& rhs) const
{ bit_field b = *this;  return b |= rhs; }

inline bit_field& bit_field::operator|= (const bit_field& rhs)
{ if(&rhs == this) return *this;
  if(size() > rhs.size()) finish = begin() + rhs.size();
  int_or(mbegin(), mend(), rhs.mbegin());
  return *this;
}

inline bit_field bit_field::operator^ (const bit_field& rhs) const
{ bit_field b = *this;  return b ^= rhs; }

inline bit_field& bit_field::operator^= (const bit_field& rhs)
{ if(&rhs == this) return *this;
  if(size() > rhs.size()) finish = begin() + rhs.size();
  int_xor(mbegin(), mend(), rhs.mbegin());
  return *this;
}

inline bit_field bit_field::operator<< (size_type t) const
{ bit_field b = *this; return b.lshift(t); }

inline bit_field& bit_field::operator<<= (size_type t)
{ return lshift(t); }

inline bit_field& bit_field::lshift(size_type t)
{
  // Note: Memcpy is 'unreliable' when doing copies of overlapping memory
  //       blocks and is therefore *not* used here.
  mem_type* b1 = mend() - (t / __WORD_BIT) - 1, * b2 = mend();
  size_type t2 = t % __WORD_BIT;
  for( ; b1 != mbegin() - 1; --b2, --b1)
  {
    if(b2 != mend()) *b2 += (*b1 >> (__WORD_BIT - t2));;
    *(b2 - 1) = (*b1 << t2);
  }
  for(--b2 ; b2 > mbegin() - 1; --b2)
    *b2 = 0;
  return *this;
}

inline bit_field bit_field::operator>> (size_type t) const
{ bit_field b = *this; return b.rshift(t); }

inline bit_field& bit_field::operator>>= (size_type t)
{ return rshift(t); }

inline bit_field& bit_field::rshift(size_type t)
{
  if(!t) return *this;

  set_zero();
  mem_type* b1 = mbegin() + (t / __WORD_BIT), * b2 = mbegin();
  size_type t2 = t % __WORD_BIT;
  for( ; b1 != mend(); ++b2, ++b1)
  {
    if(b2 != mbegin()) *(b2 - 1) += (*b1 << (__WORD_BIT - t2));
    *b2 = (*b1 >> t2); 
  }
  for( ; b2 < mend(); ++b2)
    *b2 = 0;
  return *this;
}

inline bit_field bit_field::operator~ () const
{ bit_field f = *this; return f.flip(); }

inline bit_field& bit_field::flip()
{ for(mem_type* p = mbegin(); p != mend(); ++p) (*p) = ~(*p); return *this; }

inline bit_field& bit_field::clear()
{ for(mem_type* p = mbegin(); p != mend(); ++p) (*p) = 0; return *this; }

inline bit_field& bit_field::And(const bit_field& rhs)
{
  if(&rhs == this) return *this;
  if(size() <= rhs.size()) return *this &= rhs;
  rhs.set_one();
  int_and(mbegin(), mbegin() + (rhs.mend() - rhs.mbegin()), rhs.mbegin());
  return *this;
}

inline bit_field& bit_field::Or(const bit_field& rhs)
{
  if(&rhs == this) return *this;
  if(size() <= rhs.size()) return *this |= rhs;
  rhs.set_zero();
  int_or(mbegin(), mbegin() + (rhs.mend() - rhs.mbegin()), rhs.mbegin());
  return *this;
}

inline bit_field& bit_field::Xor(const bit_field& rhs)
{
  if(&rhs == this) return *this;
  if(size() <= rhs.size()) return *this ^= rhs;
  rhs.set_zero();
  int_xor(mbegin(), mbegin() + (rhs.mend() - rhs.mbegin()), rhs.mbegin());
  return *this;
}

inline void bit_field::or_chunk(__U32 u, unsigned int offset)
{
  iterator i = start + offset;

  if(i.p > mfinish()) return;

  if(i.offset == 0)
  {
    *i.p |= u;
    return;
  }

  // Finish implementation and add similar functions elsewhere!
}

inline bit_field& bit_field::Union       (const bit_field& rhs) { return *this |= rhs; }
inline bit_field& bit_field::Intersection(const bit_field& rhs) { return *this &= rhs; }
inline bit_field& bit_field::Difference  (const bit_field& rhs) { return *this &= ~rhs; }

inline bool bit_field::subset    (const bit_field& rhs) const
{ // This does (*this & ~rhs).empty(), but more efficiently
  const_iterator final = finish;
  if(size() > rhs.size()) final = begin() + rhs.size();
  mem_type* r = rhs.mbegin();
  for(mem_type* l = mbegin(); l != final.p; ++l, ++r)
    if((*l) & ~(*r)) return 0;
  if(final.offset && ((*final.p) & ~(*r)))
  {
    const_iterator i(final.p, 0), j(r, 0);
    for(; i != final; ++i, ++j)
      if((*i) && !(*j)) return 0;
  }
  return 1;
}

inline bool bit_field::disjoint  (const bit_field& rhs) const
{ return (*this &  rhs).empty(); }

inline bool bit_field::complement(const bit_field& rhs) const
{ return (~(*this | rhs)).empty(); }

inline bool bit_field::empty() const
{
  set_zero();
  for(mem_type* p = mbegin(); p <= mfinish(); ++p)
    if(*p != 0) return 0;

  return 1;
}

inline bool bit_field::any() const
{ return !empty(); }

inline void bit_field::reserve(size_type t)
{
  if(t <= (unsigned) (mend() - mbegin()) * __WORD_BIT) return;
  bit_field b(t);

  b = *this;
  swap(b);
}

inline unsigned int* bit_field::bit_alloc(size_type n)
{
  if(n)
  {
    mem_type* m = new mem_type[(n+__WORD_BIT - 1)/__WORD_BIT];
    /*cout << "Allocating " << (n+__WORD_BIT - 1)/__WORD_BIT*4 << " bytes at "
         << m << endl;*/
    return m;
  }
  return NULL;
}

inline void bit_field::deallocate()
{ /*cout << "deallocating (dealloc): " << this << ' ' << start.p << endl;*/
  if(mbegin())
  {
    delete[] mbegin();
    start.p = NULL;
  }
}

inline void bit_field::initialize(size_type n)
{
  unsigned int* q = bit_alloc(n);
  end_of_storage = q + (n + __WORD_BIT - 1)/__WORD_BIT;
  start = iterator(q, 0);
  finish = start + n;
}

inline void bit_field::int_plus(mem_type* b1, mem_type* e1, mem_type* b2)
{
  mem_type b = 0, bt;
  for(; b1 != e1; ++b1, ++b2) { bt = *b1; *b1 += *b2 + b; b = (*b1 < bt); }
}

inline void bit_field::int_minus(mem_type* b1, mem_type* e1, mem_type* b2)
{
  mem_type b = 0, bt;
  for(; b1 != e1; ++b1, ++b2) { bt = *b1; *b1 -= (*b2 + b); b = (*b1 > bt); }
}

inline void bit_field::int_and(mem_type* b1, mem_type* e1, mem_type* b2)
{ for(; b1 != e1; ++b1, ++b2) (*b1) &= (*b2); }

inline void bit_field::int_or(mem_type* b1, mem_type* e1, mem_type* b2)
{ for(; b1 != e1; ++b1, ++b2) (*b1) |= (*b2); }

inline void bit_field::int_xor(mem_type* b1, mem_type* e1, mem_type* b2)
{ for(; b1 != e1; ++b1, ++b2) (*b1) ^= (*b2); }

inline bit_field::mem_type* bit_field::mbegin() const
{ return start.p; }

inline bit_field::mem_type* bit_field::mend() const
{ return finish.p + (finish.offset > 0); }

inline bit_field::mem_type* bit_field::mfinish() const
{ return finish.p; }

inline bit_field::size_type bit_field::bytesize() const
{ return sizeof(mem_type) * (mend() - mbegin()); }

inline bit_field::mem_type bit_field::get_mask() const
{
  mem_type mask = (1 << finish.offset) - 1; 
  return mask;
}

inline void bit_field::set_zero() const
{
  if(mfinish() && mfinish() !=mend()) *mfinish() &= get_mask();
} 
inline void bit_field::set_one()  const
{
  if(mfinish() && mfinish() != mend()) *mfinish() |= ~get_mask();
}

#endif
