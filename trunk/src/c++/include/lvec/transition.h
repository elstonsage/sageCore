#ifndef __TRANSITION_H
#define __TRANSITION_H

//
//  Recombination Transition Matrices
//
//  Author: Geoff Wedig
//
//  History:  0.1 Initial Implementation      May 22, 1998
//
//  Copyright (c) 1998 R. C. Elston
//

#include <vector>
#include "mlocus/imodel.h"
#include "lvec/meiosis_map.h"

/// \file 

/// transition.h includes two major components of the likelihood vector
/// library, the recombination_map and the marker_recombination.  These classes are
/// used with the lvector to perform Lander-Green style recombinations between
/// markers.

/// The marker_recombination is a view of the recombination_map between
/// a single pair of markers.  
class marker_recombination;

/// The recombination_map stores all the recombination fractions between
/// each pair of markers in a chromosomal region.  Each individual meiosis
/// may have a different recombination between a given set of markers, but
/// in general there is only one (or two, in the case of sex-specific
/// recombinations) To get the recombination between two markers, it is only
/// necessary to specify the first marker (the lesser one) Thus, if there
/// are n markers, there are n-1 sets of recombination fractions.
class recombination_map
{
public:

  friend class marker_recombination;

  typedef SAGE::meiosis_map    mmap;
//  typedef MultiLocus_Strategy  ml_strat;
//  typedef mmap::ind_id         ind_id;
//  typedef mmap::ped_id         ped_id;
//  typedef mmap::index          index;
//  typedef index                meiosis_type;
//  typedef mmap::size_type      size_type;
//  typedef size_type            marker_type;
//  typedef size_type            member_type;
  typedef double               recomb_frac;

  recombination_map(const meiosis_map&, );
  
  ~recombination_map() { }

// Setting values

  // Set all the values of a marker to the recombination fraction d
  r_frac set_marker  (marker_type m, r_frac d);

  // Set the meiosis r at marker m to the r. frac. d
  r_frac set_meiosis(marker_type m, meiosis_type r, r_frac d);

  // Set the mother's or father's meiosis at marker m to d
  r_frac set_mmeiosis(marker_type m, member_type i, r_frac d);
  r_frac set_fmeiosis(marker_type m, member_type i, r_frac d);
  r_frac set_mmeiosis(marker_type m, ind_id i, r_frac d);
  r_frac set_fmeiosis(marker_type m, ind_id i, r_frac d);

// Getting a value

  // Get the meiosis r's r. frac.
  r_frac meiosis(marker_type m, meiosis_type r) const;

  // Get r. frac. for mother or father of specific individual i
  r_frac mmeiosis(marker_type m, member_type i) const;
  r_frac fmeiosis(marker_type m, member_type i) const;
  r_frac mmeiosis(marker_type m, ind_id i) const;
  r_frac fmeiosis(marker_type m, ind_id i) const;

// Single Marker to Marker View object

  marker_recombination marker(marker_type m);

// Get the map this is all based on

  const mmap* meiosis_map() const;

// Size operations

  size_type  marker_count() const;
  size_type  member_count() const;
  size_type meiosis_count() const;

protected:

  typedef vector<r_frac>              r_frac_vect;
  typedef r_frac_vect::iterator       iterator;
  typedef r_frac_vect::const_iterator const_iterator;
  
  iterator       find(marker_type m, meiosis_type r);
  const_iterator find(marker_type m, meiosis_type r) const;

  LSF_ptr<mmap>     mm;
  LSF_ptr<ml_strat> mls;
  
  r_frac_vect r_fracts;

  size_type mrksize, memsize, meisize;

};

class marker_recombination
{
public:

  friend class recombination_map;

  typedef recombination_map    rmap;
  typedef rmap::mmap           mmap;
  typedef rmap::ind_id         ind_id;
  typedef rmap::meiosis_type   meiosis_type;
  typedef rmap::marker_type    marker_type;
  typedef rmap::member_type    member_type;
  typedef rmap::size_type      size_type;
  typedef rmap::r_frac         r_frac;
  
  marker_recombination() : mm(NULL) { }
  
  marker_recombination& operator= (const marker_recombination&);

// Setting values

  // Set all the values to the recombination fraction d
  r_frac set_marker  (r_frac d);

  // Set the meiosis r to the r. frac. d
  r_frac set_meiosis(meiosis_type r, r_frac d);

  // Set the mother's or father's meiosis to d
  r_frac set_mmeiosis(member_type i, r_frac d);
  r_frac set_fmeiosis(member_type i, r_frac d);
  r_frac set_mmeiosis(ind_id i, r_frac d);
  r_frac set_fmeiosis(ind_id i, r_frac d);

// Getting a value

  // Get the meiosis r's r. frac.
  r_frac meiosis(meiosis_type r) const;

  // Get r. frac. for mother or father of specific individual i
  r_frac mmeiosis(member_type i) const;
  r_frac fmeiosis(member_type i) const;
  r_frac mmeiosis(ind_id i) const;
  r_frac fmeiosis(ind_id i) const;

// Get the map this is based on

  const mmap* meiosis_map() const;

// Size operators

  size_type  member_count() const;
  size_type meiosis_count() const;

protected:

  typedef rmap::iterator iterator;
  
  marker_recombination(mmap*, iterator, iterator);

  iterator find(meiosis_type r) const;

  mmap*    mm;
  iterator f_meiosis, nf_meiosis;

};

// ========================
//     Inline Functions
// ========================

inline recombination_map::r_frac
    recombination_map::set_meiosis(marker_type m, meiosis_type r, r_frac d)
{ return *find(m, r) = d; }

inline recombination_map::r_frac
    recombination_map::set_mmeiosis(marker_type m, member_type i, r_frac d)
{ return *find(m, mm->mmeiosis(i)) = d; }

inline recombination_map::r_frac
    recombination_map::set_fmeiosis(marker_type m, member_type i, r_frac d)
{ return *find(m, mm->fmeiosis(i)) = d; }

inline recombination_map::r_frac
    recombination_map::set_mmeiosis(marker_type m, ind_id i, r_frac d)
{ return *find(m, mm->mmeiosis(i)) = d; }

inline recombination_map::r_frac
    recombination_map::set_fmeiosis(marker_type m, ind_id i, r_frac d)
{ return *find(m, mm->fmeiosis(i)) = d; }

inline recombination_map::r_frac
    recombination_map::meiosis(marker_type m, meiosis_type r) const
{ return *find(m, r); }

inline recombination_map::r_frac
    recombination_map::mmeiosis(marker_type m, member_type i) const
{ return *find(m, mm->mmeiosis(i)); }

inline recombination_map::r_frac
    recombination_map::fmeiosis(marker_type m, member_type i) const
{ return *find(m, mm->mmeiosis(i)); }

inline recombination_map::r_frac
    recombination_map::mmeiosis(marker_type m, ind_id i) const
{ return *find(m, mm->mmeiosis(i)); }

inline recombination_map::r_frac
    recombination_map::fmeiosis(marker_type m, ind_id i) const
{ return *find(m, mm->mmeiosis(i)); }

inline const SAGE::meiosis_map* recombination_map::meiosis_map() const
{ return mm; }

inline recombination_map::size_type recombination_map::marker_count() const
{ return mrksize; }

inline recombination_map::size_type recombination_map::member_count() const
{ return memsize; }

inline recombination_map::size_type recombination_map::meiosis_count() const
{ return meisize; }

inline recombination_map::iterator
    recombination_map::find(marker_type m, meiosis_type r)
{
  if(r < mm->founder_count()) return r_fracts.begin() + m * mm->founder_count() + r;
  
  return r_fracts.begin() + (marker_count() - 1) * mm->founder_count()
                          + m * mm->nonfounder_meiosis_count() + r - mmap::meiosis_bits;
}

inline recombination_map::const_iterator
    recombination_map::find(marker_type m, meiosis_type r) const
{
  if(r < mm->founder_count()) return r_fracts.begin() + m * mm->founder_count() + r;
  
  return r_fracts.begin() + (marker_count() - 1) * mm->founder_count() 
                          + m * mm->nonfounder_meiosis_count() + r - mmap::meiosis_bits;
}

inline recombination_map::marker_recombination
    recombination_map::marker(marker_type m)
{
  return marker_recombination
     (mm, find(m, 0), find(m, mmap::meiosis_bits) - mmap::meiosis_bits);
}

// ====================
// marker_recombination 
// ====================

inline marker_recombination& marker_recombination::operator=
    (const marker_recombination& mr)
{
  mm         = mr.mm;
  f_meiosis  = mr.f_meiosis;
  nf_meiosis = mr.nf_meiosis;
  
  return *this;
}

inline marker_recombination::r_frac marker_recombination::set_marker(r_frac d)
{
  for(meiosis_type r = 0; r < mm->founder_count(); ++r)
    f_meiosis[r] = d;
    
  for(meiosis_type r = 64; r < mm->nonfounder_meiosis_count() + mmap::meiosis_bits; ++r)
    nf_meiosis[r] = d;

  return d;
}

inline marker_recombination::r_frac
    marker_recombination::set_meiosis(meiosis_type r, r_frac d)
{ return *find(r) = d; }

inline marker_recombination::r_frac
    marker_recombination::set_mmeiosis(member_type i, r_frac d)
{ return *find(mm->mmeiosis(i)) = d; }

inline marker_recombination::r_frac
    marker_recombination::set_fmeiosis(member_type i, r_frac d)
{ return *find(mm->fmeiosis(i)) = d; }

inline marker_recombination::r_frac
    marker_recombination::set_mmeiosis(ind_id i, r_frac d)
{ return *find(mm->mmeiosis(i)) = d; }

inline marker_recombination::r_frac
    marker_recombination::set_fmeiosis(ind_id i, r_frac d)
{ return *find(mm->fmeiosis(i)) = d; }

inline marker_recombination::r_frac
    marker_recombination::meiosis(meiosis_type r) const
{ return *find(r); }

inline marker_recombination::r_frac
    marker_recombination::mmeiosis(member_type i) const
{ return *find(mm->mmeiosis(i)); }

inline marker_recombination::r_frac
    marker_recombination::fmeiosis(member_type i) const
{ return *find(mm->fmeiosis(i)); }

inline marker_recombination::r_frac
    marker_recombination::mmeiosis(ind_id i) const
{ return *find(mm->mmeiosis(i)); }

inline marker_recombination::r_frac
    marker_recombination::fmeiosis(ind_id i) const
{ return *find(mm->fmeiosis(i)); }

inline const SAGE::meiosis_map* marker_recombination::meiosis_map() const
{ return mm; }

inline marker_recombination::size_type marker_recombination::member_count() const
{ return mm->individual_count(); }

inline marker_recombination::size_type marker_recombination::meiosis_count() const
{ return mm->meiosis_count(); }

inline marker_recombination::marker_recombination(mmap* m, iterator i1, iterator i2)
    : mm(m), f_meiosis(i1), nf_meiosis(i2) { }

inline marker_recombination::iterator
    marker_recombination::find(meiosis_type r) const
{ if(r < mm->nonfounder_meiosis_count()) return f_meiosis + r; return nf_meiosis + r; }

#endif
