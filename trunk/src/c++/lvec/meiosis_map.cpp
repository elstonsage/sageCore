#include <cassert>
#include <algorithm>
#include "lvec/meiosis_map.h"

namespace SAGE
{

meiosis_map::meiosis_map(bool x)
           : my_filtered_mped(NULL),
             my_filtered_subped(NULL),
             my_founder_count(0),
             my_nonfounder_count(0),
             my_built(false),
             my_x_linked(x)
{}

meiosis_map::meiosis_map(const FPED::Multipedigree& fmp, bool x)
           : my_filtered_mped(&fmp),
             my_filtered_subped(NULL),
             my_founder_count(0),
             my_nonfounder_count(0),
             my_built(false),
             my_x_linked(x)
{}

meiosis_map::meiosis_map(subpedigree_const_pointer sp, bool x)
           : my_filtered_mped((multipedigree*) sp->multipedigree()),
             my_filtered_subped(sp),
             my_founder_count(0),
             my_nonfounder_count(0),
             my_built(false),
             my_x_linked(x)
{
  if( sp )
  {
    set_subpedigree(sp);
    build();
  }

  if( built() ) set_masks();
}

meiosis_map::meiosis_map(const meiosis_map& mm)
           : my_filtered_mped(mm.my_filtered_mped),
             my_filtered_subped(mm.my_filtered_subped),
             my_founder_count(mm.my_founder_count),
             my_nonfounder_count(mm.my_nonfounder_count),
             my_built(mm.my_built),
             my_x_linked(mm.my_x_linked)
{
  if( built() ) set_masks();
}

meiosis_map&
meiosis_map::operator=(const meiosis_map& mm)
{
  my_filtered_mped    = mm.my_filtered_mped;
  my_filtered_subped  = mm.my_filtered_subped;
  my_founder_count    = mm.my_founder_count;
  my_nonfounder_count = mm.my_nonfounder_count;

  my_built  = mm.my_built;
  my_x_linked = mm.my_x_linked;

  p_meioses = mm.p_meioses;
  masks     = mm.masks;
  fm        = mm.fm;
  nfm       = mm.nfm;
  fbits     = mm.fbits;
  nfbits    = mm.nfbits;

  return *this;
}

bool
meiosis_map::unbuild()
{
  if( !my_built ) return true;

  masks.resize(0);

  my_built = false;

  return true;
}

meiosis_map::subpedigree_const_pointer 
meiosis_map::set_subpedigree(subpedigree_const_pointer sp)
{
  masks.resize(0);

  my_filtered_subped = sp;

  my_built = false;

  return my_filtered_subped;
}

bool
meiosis_map::build()
{
  if( built() ) return true;

  my_founder_count    = 0;
  my_nonfounder_count = 0;
  
  multipedigree::member_const_iterator ind = my_filtered_subped->member_begin();

  for( ; ind != my_filtered_subped->member_end(); ++ind )
  {
    if( (!ind->parent1() || !ind->parent2()) && ind->offspring_count() )
      ++my_founder_count;

    if( ind->parent1() && ind->parent2() )
      ++my_nonfounder_count;
  }

  return my_built = true;
}

bool
meiosis_map::set_masks()
{
  fbits = 0;

  p_meioses.resize(0);
  p_meioses.resize(my_filtered_subped->member_count());

  // We assign a number for each nonfounder bit such that each founder bit
  // gets a proper value
  nfbits = meiosis_bits;

#if 0
  cout << endl << "meiosis_map::set_masks()... " << endl
       << "start : "
       << "fbits = " << fbits << ", nfbits = " << nfbits << ", meiosis_bits = " << meiosis_bits
       << ", fm = " << fm << ", nfm = " << nfm << endl;
#endif

  for( size_type i = 0; i < my_filtered_subped->member_count(); ++i )
  {
    if( founder(i) )         // If founder
    {
      masks.push_back(0);
    }
  }

  for( size_type i = 0; i < my_filtered_subped->member_count(); ++i )
  {
#if 0
  cout << "i = " << i << ", " << my_filtered_subped->member_index(i).name() << endl;
#endif

    if( !founder(i) )         // If founder
    {
      index mom = mother_index(i);
      index dad = father_index(i);

#if 0
  cout << "mom = " << mom << ", dad = " << dad << ", "
       << p_meioses[i].moth << ", " << p_meioses[i].fath << endl;
#endif

      set_index(mom, p_meioses[i].moth);
      set_index(dad, p_meioses[i].fath);
    }
  }

  // Since we used the mother of each founder to keep track of the mask
  // index for that founder, it is set to the first child's index.  We
  // must reset this to -1 now.

  for( size_type i = 0; i < my_filtered_subped->member_count(); ++i )
  {
    if( founder(i) )
      p_meioses[i].moth = index_err;
  }

#if 0
  cout << "middle : "
       << "fbits = " << fbits << ", nfbits = " << nfbits << ", meiosis_bits = " << meiosis_bits
       << ", fm = " << fm << ", nfm = " << nfm << endl;
#endif

  nfbits -= meiosis_bits;

  // Now we create our masks and we're done.

  fm  = (((storage_type) 1) <<    founder_meiosis_count()) - 1;
  nfm = (((storage_type) 1) << nonfounder_meiosis_count()) - 1;

#if 0
  cout << "end : "
       << "fbits = " << fbits << ", nfbits = " << nfbits << ", meiosis_bits = " << meiosis_bits
       << ", fm = " << fm << ", nfm = " << nfm << endl;

  dump_map();
#endif

  return true;
}

void
meiosis_map::set_index(size_type p, index& c)
{
#if 0
  cout << "set_index()..with p index = " << p << ", p meiosis = " << c << endl;
#endif

  if(p_meioses[p].moth == index_err)           // if first child
  {
    // We use moth of the parent (p) as a reference into the masks vector
    // while doing our construction.  It will be reset to -1 when the
    // construction is complete

    p_meioses[p].moth = c = fbits++;   // Set index

#if 0
  cout << "p_meioses[" << p << "].moth = " << p_meioses[p].moth << ", fbits = " << fbits << endl;
#endif

  }
  else                                         // Otherwise
  {
    c = nfbits++;

#if 0
  cout << c << ", nfbits = " << nfbits << endl;
#endif

    // If the parent is a founder, we add this child to the mask.
    if( founder(p) )
#if 0
{
  cout << "masks[ p_meioses[" << p << "].moth ] = " << masks[ p_meioses[p].moth ] << endl;
#endif

      masks[ p_meioses[p].moth ] ^= ((storage_type) 1) << (c - meiosis_bits);

#if 0
  cout << "p_meioses[" << p << "].moth = " << p_meioses[p].moth
     << ", masks[ p_meioses[" << p << "].moth ] = " << masks[ p_meioses[p].moth ] << endl;
}
#endif
  }
}

void
meiosis_map::dump_map() const
{
  cout << "p_meioses :" << endl;
  for( size_t i = 0; i < p_meioses.size(); ++i )
  {
    cout << i << " " << member(i)->name() <<  " : "
         << p_meioses[i].moth << ", " << p_meioses[i].fath
         << ", " << mask(p_meioses[i].moth) << ", " << mask(p_meioses[i].fath)
         << endl;
  }

  cout << endl;
  cout << "masks :" << endl;
  for( size_t i = 0; i < masks.size(); ++i )
  {
    cout << i << " : " << masks[i] << endl;
  }

  cout << endl;
  cout << "fm = "    << fm    << ", nfm = " << nfm << endl;
  cout << "fbits = " << fbits << ", nfbits = " << nfbits << endl;
  cout << "founder_count = " << my_founder_count
       << ", nonfounder_count = " << my_nonfounder_count << endl;

  cout << (storage_type) 1;
}

}
