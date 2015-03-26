#ifndef __INHERITANCE_VECTOR_H
#define __INHERITANCE_VECTOR_H

//
//  Inheritance Vectors and Descent Trees - Descent based information
//     structures.
//
//  Author: Geoff Wedig
//
//  History:  0.1 Initial Implementation      Sept 25 1997
//            0.2 Revised, added to, and      Mar  24 1998
//                restructured
//            1.0 Updated to new library      Jul     2001
//
//  Copyright (c) 2001 R. C. Elston
//

#include <vector>
#include "LSF/LSF.h"
#include "lvec/meiosis_map.h"

namespace SAGE
{

class iv_reference
{
public:

  friend class inheritance_vector;

  typedef meiosis_map::index        index;
  typedef meiosis_map::mask_type    mask_type;
  typedef meiosis_map::storage_type storage_type;
  typedef iv_reference              reference;

  iv_reference() : p(NULL), eq(NULL), mask(0), eqmask(0) { }

// Conversion and copy operators.

  inline operator bool() const;

  inline reference& operator=(bool x);
  inline reference& operator=(const reference& x);

// Boolean functions

  inline bool operator== (const reference& x) const;
  inline bool operator<  (const reference& x) const;

  inline void flip();

protected:

  iv_reference(storage_type* d, storage_type* q, mask_type m, mask_type n)
      : p(d), eq(q), mask(m), eqmask(n) { }

  storage_type* p, * eq;
  mask_type     mask, eqmask;

};

/// Stores the bits describing the descent of alleles.
/** The inheritance_vector shares many features with integers and iterators.
 *  Its main job is to store the bits describing the descent of alleles
 *  through each meiosis in a pedigree.  It supports a number of the iterator
 *  functions (++, --, etc) as well as features allowing setting and
 *  modification of bits.
 * 
 *  There are two methods of iteration through the inheritance_vector values.  In
 *  the first (default) we iterate through the values of the equivalence
 *  class.  In the second, we iterate through all the values the
 *  inheritance_vector can take.  In both, the isEnd() function is true if
 *  we pass the end of the valid values of the inheritance_vector.
 */
class inheritance_vector
{
public:

  /// Iteration method
  enum iter_mthd
  {
    EQ,        //< Equivalence Class
    BW         //< Bitwise
  };

  typedef meiosis_map::storage_type         storage_type;
  typedef meiosis_map::mask_type            mask_type;
  typedef meiosis_map::index                index;
  typedef meiosis_map::size_type            size_type;
  typedef iv_reference                      reference;
  typedef bool                              const_reference;
  typedef storage_type                      equivalence_class;
  typedef meiosis_map::member_const_pointer member_const_pointer;

  inheritance_vector(iter_mthd i = EQ)
      : mm(), found(0), nonf(0), eq(0), iter(i) { }
  explicit inheritance_vector(const meiosis_map& m, iter_mthd i = EQ)
      : mm(), iter(i)
  { (*this)(m); }

  ~inheritance_vector() { }
  
  inheritance_vector& operator=  (const inheritance_vector& i);
  inheritance_vector& operator() (const meiosis_map& m);
  
  const meiosis_map& get_meiosis_map() const { return mm; }

  inline storage_type      get_founders()          const;
  inline storage_type      get_nonfounders()       const;
  inline equivalence_class get_equivalence_class() const;
  
  inline equivalence_class set_equivalence_class   (equivalence_class s);
  inline size_type         num_equivalence_classes () const;

  /// @name Getting individual bits
  //@{
  inline reference       meiosis(index i);
  inline const_reference meiosis(index i) const;

  inline reference       operator[] (index i);
  inline const_reference operator[] (index i) const;

  inline reference       mother_meiosis(member_const_pointer i);
  inline const_reference mother_meiosis(member_const_pointer i) const;
  inline reference       father_meiosis(member_const_pointer i);
  inline const_reference father_meiosis(member_const_pointer i) const;
  //@}

  /// @name Operator ++ --, etc
  //@{
  inline iter_mthd get_iteration_method() const;
  inline iter_mthd set_iteration_method(iter_mthd i);

  inheritance_vector& operator++ ();
  inheritance_vector& operator-- ();

  inheritance_vector  operator++ (int);
  inheritance_vector  operator-- (int);
  //@}

  bool isBegin() const;
  bool   isEnd() const;

protected:

  storage_type inc_founder();
  storage_type dec_founder();

  meiosis_map        mm;                //< My meiosis map

  storage_type       found, nonf;       //< storage of inheritance bits
  equivalence_class  eq;                //< Equivalence class storage

  iter_mthd          iter;              //< Iteration method
};


#include "lvec/inheritance_vector.ipp"

} // End namespace

#endif
