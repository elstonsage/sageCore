#ifndef __FFT_BIT_COUNT_H
#define __FFT_BIT_COUNT_H

#include "lvec/meiosis_map.h"
#include "lvec/inheritance_vector.h"

namespace SAGE
{

class fft_bit_count
{
public:

  typedef inheritance_vector::equivalence_class equivalence_class;

// Constructors/Destructor
  fft_bit_count();

  ~fft_bit_count();
  
  bool build(long bits);
  
  bool count(meiosis_map*);
  
// Accessors

  meiosis_map* get_meiosis_map() const;

  long operator[](equivalence_class) const;
  
  equivalence_class capacity() const;
  equivalence_class member_count() const;

protected:

  void deallocate();

  void initialize(equivalence_class n);
  
// data members

  vector<int> storage;

  bool built;
  
  meiosis_map* mm;
};

#include "lvec/fft_bit_count.ipp"

}

#endif
