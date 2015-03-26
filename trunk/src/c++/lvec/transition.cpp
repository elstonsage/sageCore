#include "lvec/transition.h"

recombination_map::recombination_map(mmap* _mm, ml_strat* _mls)
   : mm(_mm), mls(_mls)
{
  // Check to make sure everything is ok.
  if(!mm || !mls || !mm->good() || !mls->good())
  {
    setstate(badbit);

    mm  = NULL;
    mls = NULL;

    return;
  }

  mrksize = mls->size();
  memsize = mm->individual_count();
  meisize = mm->meiosis_count();

  // Number of recombination fractions is the # of meioses * (# of markers - 1)
  r_fracts.resize(meiosis_count() * (marker_count() - 1));

  for(marker_type m = 0; m < marker_count()-1; ++m)
    set_marker(m, mls->dist(m));
}

/// Set all the values of recombination fraction between marker m and marker
/// m+1 to the recombination fraction d

recombination_map::r_frac recombination_map::set_marker
    (marker_type m, r_frac d)
{
  iterator i = find(m, 0);

  for(meiosis_type r = 0; r < mm->founder_count(); ++r) i[r] = d;

  i = find(m, mmap::meiosis_bits);

  for(meiosis_type r = 0; r < mm->nonfounder_meiosis_count(); ++r) i[r] = d;

  return d;
}
