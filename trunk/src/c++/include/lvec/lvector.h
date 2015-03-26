#ifndef __LVECTOR_H
#define __LVECTOR_H

//
//
//  Likelihood Vector 0.0 - Storage of inheritance vector likelihoods
//
//  Author:  Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History: 0.1 gcw Revised based on bit field            May  5 1998
//           1.0 gcw Final version for initial 4.0 Release Dec    2001
//           1.1 gcw Upgraded to new libraries             May    2002
//
//  Note: Do not assume that vector operations work.
//  
//  Copyright (c) 1998-2002  R.C. Elston

//#include <utility>
//#include <algobase>

#include "numerics/log_double.h"
#include "mlocus/imodel.h"
#include "lvec/inheritance_vector.h"
// Removed this capability for now.  We don't need it yet.
//#include "lvec/transition.h"
#include "lvec/fft_bit_count.h"
#include "lvec/lvec_allocator.h"

#include <iomanip>

namespace SAGE
{

class const_lvector_iterator;
class Likelihood_Vector;

class const_lvector_iterator
{
public:

  friend class Likelihood_Vector;

  typedef Likelihood_Vector                     lvector;
  typedef const_lvector_iterator                const_iterator;
  typedef double                                likelihood;
  typedef inheritance_vector::equivalence_class equivalence_class;
  typedef likelihood                            const_reference;
  typedef likelihood                            value_type;
  typedef equivalence_class                     difference_type;
  typedef size_t                                size_type;

// Constructors

  const_lvector_iterator() : p(0) {}
  const_lvector_iterator(const const_iterator& b) : p(b.p) {}

// Dereferencing

  const_reference operator*() const;

  const_reference operator[](difference_type i) const;

// Increment/Decrement

  const_iterator& operator++();
  const_iterator  operator++(int);
  const_iterator& operator--();
  const_iterator  operator--(int); 

  const_iterator& operator+=(difference_type i);
  const_iterator& operator-=(difference_type i);

  const_iterator operator+(difference_type i) const;
  const_iterator operator-(difference_type i) const;

  difference_type operator-(const_iterator x) const;

// Comparisons

  bool operator==(const const_iterator& x) const;
  bool operator!=(const const_iterator& x) const;

  bool operator<(const const_iterator& x) const;

protected:

  const_lvector_iterator(const likelihood* x) : p(x) {}

  const likelihood* p;
};

class Likelihood_Vector
{
public:

  friend class lv_acceptor;

  typedef Likelihood_Vector                     lvector;
  typedef double                                likelihood;
  typedef inheritance_vector::equivalence_class equivalence_class;
  typedef likelihood                            const_reference;
  typedef likelihood                            value_type;
  typedef equivalence_class                     difference_type;
  typedef size_t                                size_type;  
  typedef const_lvector_iterator                const_iterator;
  typedef meiosis_map                           mmap;
  typedef fft_bit_count                         fft;

// Constructors/Destructor
  Likelihood_Vector();
  explicit Likelihood_Vector(size_type n);
  Likelihood_Vector(mmap*, const SAGE::MLOCUS::inheritance_model&);
  Likelihood_Vector(const lvector&);
//  Likelihood_Vector(const lvector&, const Marker_Recombination&);
  Likelihood_Vector(const lvector&, const mmap*, double theta);

  ~Likelihood_Vector();

  // Pseudo-constructor.  Generates the lv with the mmap and pm given
  lvector& operator() (mmap*, const SAGE::MLOCUS::inheritance_model&, bool use_pop_freq = true);

// operator =

  lvector& operator=  (const lvector& rhs);

// Memory Management

  // Get rid of all extra memory (ie, if maxbits > bits)
  // Currently not available.
//  void compress();

// Multiplication style operations

  // Multiply by the Marker_Recombination
//  lvector& operator*= (const Marker_Recombination&);

  // Multiply by the transition matrix formed by the mmap and r.f. theta
  lvector& operator() (const mmap*, double theta);

  // Multiply using fft method using the fft_bit_count and r.f. theta
  lvector& operator() (const fft&, double theta);

  // This multiplies each element by each element
  lvector& operator*= (const lvector&);

  // This scales the vector by the log_double
  lvector& operator*= (log_double);

  // This adds two lvs together 
  lvector& operator+= (const lvector&);

// Data modification and retrieval options

  // Flatten sets to a 1/2^n vector, each element being 1/2^n
  void flatten(long n);

  void normalize();
  
  likelihood total()     const;
  log_double log_total() const; // Returns the total * scale term

  log_double log_scale() const; // Returns the scale term

  double     information() const;

  equivalence_class fixed_bits() const { return fixed; }

// Iterators and references

  const_reference operator[] (size_type n) const;

  const_iterator begin() const;
  const_iterator   end() const;
  
  const_reference front() const;
  const_reference  back() const;
  
// Miscellaneous operations

  size_type         size() const;
  size_type     capacity() const;
  size_type    bit_count() const;
  size_type bit_capacity() const;


  void swap(lvector& rhs);
  
// Operators

  // Equality operators

  bool operator== (const lvector& rhs) const;
  bool operator!= (const lvector& rhs) const;

  bool is_valid() const;
  bool set_valid(bool);

protected:

#if defined(BROKEN_VECTOR_ALLOCATOR)
  typedef vector<likelihood> storage_type;
#else
  typedef vector<likelihood, lvec_allocator> storage_type;
#endif

  void increment_value(equivalence_class n, equivalence_class dont_care_bits, likelihood p);

  equivalence_class  start() const;
  equivalence_class& bump_up(equivalence_class&) const;
  
  void deallocate();

  void clear_bits();
  void resize(size_type n);
  void resize_and_clear_bits(size_type n);
  
// For multiplication

  void Process(const mmap*, vector<double>::iterator, double);

  void Process_Fixed   (equivalence_class, size_t, vector<double>::iterator, const mmap*);
  void Process_Natural (equivalence_class, size_t, vector<double>::iterator, const mmap*);

  void Process_Founders(const mmap*, vector<double>::iterator);

// data members

  storage_type storage;
  
  equivalence_class first, fixed;

  bool        f;
  size_type   bits, max_bits;
  log_double  my_scale;

  bool        my_valid;
};

#include "lvec/lvector.ipp"

}

#endif
