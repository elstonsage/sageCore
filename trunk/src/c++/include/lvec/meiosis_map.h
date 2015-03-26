#ifndef MEIOSIS_MAP_H
#define MEIOSIS_MAP_H

//==========================================================================
//  File:    meiosis_map.h
//
//  Author:  Geoff Wedig
//
//  History: 0.1 Initial Implementation                         Sept 25 1997
//           0.2 Revised and moved into own file                Mar  24 1998
//           0.3 Split into Marker specific                     Apr  07 1998
//               and non-marker spec. parts
//           1.0 Modified to be independent of pedigree_section yjs Dec.2003
//               and utilize RPED::filtered_multipedigree     
//           1.1 Added X-linkage                                yjs Apr.2004
//
//  Notes:   Meiosis Mapping - Ordered meioses and individuals for a
//           pedigree to generate likelihood vectors.
//
//  Copyright (c) 1998 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include <boost/integer.hpp>
#include <vector>
#include "LSF/LSF.h"
#include "fped/fped.h"

namespace SAGE
{


#define BIT_COUNT 32

/// Storage in the meiosis_map requires a fixed, understood number of bits. 
/// This is an unsigned 64 bit integer.  That is larger than is needed (32
/// would probably be adequate for the forseeable future) in some cases, but
/// is useful in that it provides something that we're unlikely to need
/// changing later and does not cost us much.

typedef boost::uint_t<BIT_COUNT>::fast meiosis_map_storage_type;

/// Extends the pedigree_section to include and order meioses

/// The meiosis_map extends the pedigree_section into a set of mappings from individual to
/// a bit in an inheritance vector.  This mapping makes it possible to
/// perform many of the tasks necessary in Likelihood_Vector generation. 
/// Due to the complexity of these operations, the meiosis_map only operates
/// on small to medium sized pedigrees.
///
/// There exist within the meiosis map symmetries which can be taken
/// advantage of.  Founder symmetries are those symmetries that exist due to
/// founder phase being unknown.  In these cases, the bit of the first child
/// descended from a particular founder is specified as a founder bit.  A
/// founder bit has several properties that reduce the size of the
/// Likelihood_Vector.
///
/// The Functions it provides above those provided by the pedigree_section are:
///
///   -#  Id, sex -> int.  This is the reference of the allele into an
///       inheritance vector, returning -1 for founders.
///
///       Indicies for first meioses of founders are on the range [0..f) and
///       indices for non-founders on the range [x..x+2n-f) (where x is 64
///       bits).
///
///   -#  int -> storage_type.  This storage_type is the type mask used in
///       the equivalence class.  The int passed is the reference of the
///       allele changing.  For non-founder bits, this is simply the bit to
///       be changed.  For founder bits (only applicable when founder
///       symmetries are used), it is the bits of the entire sibship.
///
/// These functions are provided by two internal data structures.  The first
/// two functions are supported by a data structure storing each individual
/// with two indices for the meiosis from the mother and father.  The second
/// structure is a vector of masks for each founder.  The mask value returned
/// is the element for which the int is the index into this vector (ie, the
/// founder value) if its the first child of a founder, and 1 << (index) if
/// not.

class meiosis_map
{
public:

  typedef FPED::Multipedigree                      multipedigree;

  typedef multipedigree::member_type               member_type;
  typedef multipedigree::member_pointer            member_pointer;
  typedef multipedigree::member_const_pointer      member_const_pointer;
  typedef multipedigree::pedigree_type             pedigree_type;
  typedef multipedigree::pedigree_pointer          pedigree_pointer;
  typedef multipedigree::pedigree_const_pointer    pedigree_const_pointer;
  typedef multipedigree::subpedigree_type          subpedigree_type;
  typedef multipedigree::subpedigree_pointer       subpedigree_pointer;
  typedef multipedigree::subpedigree_const_pointer subpedigree_const_pointer;

  typedef unsigned int                             index;
  typedef unsigned int                             size_type;
  typedef std::list<index>                         individual_list;
  typedef individual_list::const_iterator          individual_iterator;

  typedef meiosis_map_storage_type  storage_type;
  typedef storage_type              mask_type;

  static const index index_err    = (index) -1;
  static const index meiosis_bits = BIT_COUNT;

  /// @name Constructors/Destructor
  //@{
  meiosis_map (bool x = false);
  meiosis_map (const multipedigree& mped, bool x = false);
  meiosis_map (subpedigree_const_pointer  ped, bool x = false);
  meiosis_map (const meiosis_map&);
  
  ~meiosis_map() { }

  meiosis_map& operator=(const meiosis_map&);
  //@}

  const multipedigree&          get_multipedigree() const;
  pedigree_const_pointer        get_pedigree() const;
  subpedigree_const_pointer     get_subpedigree() const;

  /// @name Construction objects
  //@{

  /// Set up the meiosis map for changes
  bool unbuild();

  /// Change the subpedigree_const_pointer.  This calls unbuild as well.
  subpedigree_const_pointer set_subpedigree(subpedigree_const_pointer);

  /// Build builds the internal data structures.  If there are no individuals
  /// specified, it uses all individuals of the current subpedigree.
  bool build();

  // returns if there has been a build since set_ped() and if that build was
  // successful.
  bool built() const;
  //@}

  /// @name Mother and Father meiosis positions in Inheritance Vector
  //@{
  index mother_meiosis (member_const_pointer id) const;
  index father_meiosis (member_const_pointer id) const;

  index mother_meiosis (size_type i) const;
  index father_meiosis (size_type i) const;
  //@}
  
  /// @name Masks
  //@{
  const mask_type&    founder_mask() const;
  const mask_type& nonfounder_mask() const;

  mask_type mask(index i) const; 

  mask_type mother_mask(size_type i) const;
  mask_type father_mask(size_type i) const;

  mask_type mother_mask(member_const_pointer i) const;
  mask_type father_mask(member_const_pointer i) const;
  //@}
  
  // @name Number of meioses (total, founder and non-founder)
  //@{
  size_type meiosis_count           () const;
  size_type founder_meiosis_count   () const;
  size_type nonfounder_meiosis_count() const;
  //@}
  
  // Get id status
  //
  bool founder(member_const_pointer id) const;
  bool founder(size_type             i) const;

  member_const_pointer member(size_type i) const;
  member_const_pointer operator[](size_type i) const;

  member_const_pointer mother(member_const_pointer id) const;
  member_const_pointer father(member_const_pointer id) const;

  member_const_pointer mother(size_type             i) const;
  member_const_pointer father(size_type             i) const;

  size_type mother_index(member_const_pointer id) const;
  size_type father_index(member_const_pointer id) const;

  size_type mother_index(size_type             i) const;
  size_type father_index(size_type             i) const;

  // Number of founders and nonfounders
  size_type founder_count() const;
  size_type nonfounder_count() const;

  size_t    family_count() const;
  size_t    individual_count() const;

  // 2n-f where n is nonfounders and f is founders
  size_type bit_count() const;

  bool is_x_linked() const;
  bool is_father_bit(index) const;

  void dump_map() const;

protected:

  struct parent_meioses
  {
    index moth;
    index fath;
    
    parent_meioses() : moth(index_err), fath(index_err) { } 
  };

  /// Construction functions
  //@{
  bool set_masks();
  void set_index(size_type p, index& c);
  //@}
  
  std::vector<parent_meioses> p_meioses;
  std::vector<mask_type>      masks;
  mask_type                   fm, nfm; //< Founder and non-founder masks
  size_type                   fbits;   //< Number of founder meioses
  size_type                   nfbits;  //< Number of non-founder meioses

  const multipedigree*        my_filtered_mped;
  subpedigree_const_pointer   my_filtered_subped;
  size_type                   my_founder_count;
  size_type                   my_nonfounder_count;
  bool                        my_built;
  bool                        my_x_linked;
};

#undef BIT_COUNT

#include "lvec/meiosis_map.ipp"

}

#endif
