#include "lvec/fft_bit_count.h"

namespace SAGE
{

bool fft_bit_count::count(meiosis_map* _mm)
{
  if(!_mm) return false;

  if(_mm == mm) return true;

  mm = _mm;

  initialize(mm->nonfounder_meiosis_count());

  for(equivalence_class i = 1; i < storage.capacity(); ++i)
  {
    storage[i] = 0;

    for(equivalence_class j = 1; j < storage.capacity(); j <<= 1)
      if(i & j) ++storage[i];

    int bits = 0;

    for(size_t j = 0; j < mm->founder_count(); ++j)
    {
      equivalence_class e = i & mm->mask(j);
      bits = 0;

      for(size_t k = 1; k < storage.capacity(); k <<=1)
        bits += (bool) (e & k);

      if(bits % 2) ++storage[i];
    }
  }

  return true;
}

}
