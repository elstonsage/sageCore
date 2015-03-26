#ifndef MCMC_MEIOSIS_MAP_H
#define MCMC_MEIOSIS_MAP_H

//==========================================================================
//  File:    mcmc_meiosis_map.h
//
//  Author: Geoff Wedig
//          Yeunjoo Song
//
//  History:  0.1 Initial Implementation                        Sept 25 1997
//            0.2 Revised and moved into own file               Mar  24 1998
//            0.3 Split into Marker specific                    Apr  07 1998
//                and non-marker spec. parts
//
//            1.0 Updated to new libraries                       yjs May. 04
//
//  Notes:   Meiosis Mapping - Ordered meioses and individuals for a
//           pedigree to generate likelihood vectors.
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/definitions.h"

namespace SAGE
{

namespace MCMC
{

/// \brief Manages the meiosis indices for larger pedigrees, particularly for MCMC simulation
///
/// The McmcMeiosis map assigns indices for each meiosis in a single, connected
/// subpedigree.  Each non-founder in the pedigree has two meioses, one from
/// the mother, one from the father.
///
class McmcMeiosisMap
{

public:

  /// \name Object Management
  //@{
  McmcMeiosisMap (const FPED::Subpedigree&   sped,  bool x = false);
  McmcMeiosisMap (const McmcMeiosisMap&);
  
  ~McmcMeiosisMap() { }

  McmcMeiosisMap& operator=(const McmcMeiosisMap&);
  //@}

  const FPED::Multipedigree::multipedigree_type&  get_multipedigree() const;
  const FPED::Subpedigree&                        get_subpedigree()   const;

  /// \name Mother and Father meiosis positions in Inheritance Vector
  //@{
  size_t get_mother_meiosis(const FPED::Member& id) const;
  size_t get_father_meiosis(const FPED::Member& id) const;

  //@}
  
  /// \name Number of meioses (total, founder and non-founder)
  //@{
  size_t get_meiosis_count           () const;
  //@}

  // Get id status
  //
  FPED::MemberConstPointer get_mother(const FPED::Member& id) const;
  FPED::MemberConstPointer get_father(const FPED::Member& id) const;

  size_t get_founder_count()    const;
  size_t get_nonfounder_count() const;

  const FPED::Member& get_child_of_meiosis(size_t i) const;
  
  size_t get_individual_count()      const;
  size_t get_family_count()          const;

  bool is_x_linked()             const;

  void dump_map(ostream& o)      const;

protected:

  static const size_t index_err = (size_t) -1;

  /// Construction functions
  //@{
  void count_members();
  void set_indices();
  //@}
  
  /// \brief stores the indices of the mother and father meiosis bits for an individual
  struct parent_meioses
  {
    size_t moth;
    size_t fath;
    
    parent_meioses() : moth(index_err), fath(index_err) { } 
  };

  std::vector<parent_meioses> my_parent_meioses;       ///< Storage of the meioses for each non-founder
  std::vector<size_t>         my_individual_positions; ///< Given a meiosis, indexes which individual is the child

  FPED::Multipedigree::multipedigree_const_pointer my_mped;
  FPED::SubpedigreeConstPointer                    my_subped;

  size_t                      my_meiosis_count;        ///< Number of meioses

  size_t                      my_founder_count;
  size_t                      my_nonfounder_count;

  bool                        my_x_linked;
};

#include "mcmc/mcmc_meiosis_map.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif

