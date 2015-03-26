#ifndef __GENO_ELIMINATE_H
#define __GENO_ELIMINATE_H

//
// Genotype Elimination using the Pedigree_Marker (ped_marker.h)
//
// Author: Geoff Wedig (wedig@darwin.cwru.edu
//
// History: 0.1 gcw Initial Implementation
//          1.0 gcw Rewritten and Revised for RefMultiPedigree
//          2.0 yjs Added X & Y-linkage  - Mar. 2002
//
//
//
// Copyright (c) 2001 by R. C. Elston

#include <vector>
#include <utility>
#include "LSF/LSF.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "fped/fped.h"
#include "containers/bitfield.h"
#include "mlocus/imodel.h"
#include "gelim/inconsistency_handler.h"
#include "gelim/valid_parental_genotypes.h"

#ifdef __KCC
using std::pair;
using std::vector;
#endif

namespace SAGE
{

/// Genotype_Eliminator performs genotype elimination on an inheritance
/// model for a particular subpedigree. It calls an inconsistency_handler
/// object upon finding inconsistent individuals.  When perfoming checks on
/// families that have inconsistent individuals, families involving the
/// inconsistent individual as a parent are ignored while the individual is
/// ignored in the nuclear family to which it is a child.  This is to avoid
/// propogating errors throughout the subpedigree, but still allow as much
/// checking as is feasible.  Note that if a child is inconsistent with the
/// parents, the parent will be labeled inconsistent.  Children are only
/// labeled inconsistent if they are not consistent with parental genotypes,
/// not if siblings are inconsitent with the parental genotypes.  This is to
/// avoid error propogation.  As an example, in the following subpedigree,
///
/// \verbatim
///
///     1/1 --- 1/2
///         _|_
///        |   |
///       2/2 1/2
///
/// \endverbatim
///
/// both parents and the first child would be labeled inconsistent, but the
/// second would not be.
///
/// If there are errors before doing genotype elimination, genotype
/// elimination is not performed.  If there are errors in the subpedigree, the
/// Pedigree_Marker's fail bit is set to true upon return.

class genotype_eliminator
{
public:

  enum removal_type { none, all, genotype };

  typedef MLOCUS::inheritance_model         imodel;
  typedef FPED::Multipedigree               multipedigree_type;
  typedef FPED::Pedigree                    pedigree_type;
  typedef FPED::Subpedigree                 subpedigree_type;
  typedef FPED::Family                      family_type;
  typedef FPED::Member                      individual_type;
  typedef const individual_type*            individual_const_pointer;
  typedef inconsistency_handler             handler;
  typedef handler::error_type               error_type;

  genotype_eliminator();
  
  ~genotype_eliminator() { }

  // Change the inconsistency_handler element
  
  const handler& get_errors() const { return err; }

  void clear_errors() { err.clear(); }

  void mark_uninformative(imodel& model, size_t marker);

  bool set_subpedigree(const subpedigree_type&);

//  bool invalidate_individual(const individual_type&);

  /// \brief Processes the subpedigree and returns an error code.
  ///
  /// Note, the imodel here and elsewhere is assumed to be specific to the
  /// subpedigree that has been set, and that the model and the section
  /// have a 1<->1 mapping.  If this is not the case, badness may result.
  ///
  /// Return types as follows:
  /// 0 - Everything good, genotype elimination performed fine.
  /// 1 - Model uninformative.  No GE performed.
  /// 2 - Model inconsistency when doing GE.
  size_t process(imodel&, size_t marker, removal_type = genotype, bool do_x_data_check = true);

  /// \brief Does genotype elimination on a single nuclear family.
  ///        specified by a sibling in the family.
  ///
  /// removal_type and propagate control what information, if any, is
  /// transfered to surrounding nuclear families.
  ///
  /// Return types as follows:
  /// 0 - GE done, everything ok.
  /// 1 - Invalid family.  No family with those parents (may mean no
  ///     children shared by those people) in subpedigree.
  /// 2 - Parental Genotypes already invalid.  GE cannot be performed.
  /// 3 - Family inconsistent.
  size_t process_family(imodel&, size_t marker, const family_type& fam,
                        error_type e, removal_type = genotype,
                        bool propogate = true);

protected:

  typedef const individual_type*            ind_id;

  typedef list<const family_type*> fam_set;   // Set of families, indexed by first child.

  typedef imodel::phased_penetrance_iterator    phased_iterator;   // Genotype iterator
  typedef imodel::unphased_penetrance_iterator  unphased_iterator; // Genotype iterator
  typedef bit_field                             bit_type;          // The bit class (see note below)

  struct child_info_type
  {
    individual_const_pointer ptr;
    bit_type                 data;
    phased_iterator          begin;
    phased_iterator          end;
  };

  typedef vector<child_info_type>               child_vector;      // Child vector

  typedef valid_parental_genotypes                     par_genotypes;
  typedef par_genotypes::parental_genotype_pair_vector par_pair_vector;
  typedef par_genotypes::genotype_pair_iterator        par_pair_iterator;
  typedef par_genotypes::genotype_iterator             par_iterator;

  void generate_list(imodel&); // Generate the initial list of fams

  // Pass 1.
  void generate_valid_parental_genotypes   (imodel&, par_genotypes&);

  void generate_all_parental_genotypes     (imodel&, par_pair_vector&);

  // Pass 2.
  void generate_valid_child_genotypes    (imodel&, par_pair_iterator begin, par_pair_iterator end);

  // Error Detection - returns true if errors detected 
  bool generate_errors  (imodel&, const par_genotypes&, size_t, error_type);

  // Updating data
  void remove_genotypes(imodel&, const par_genotypes&, bool propagate);
  void remove_all      (imodel&, const par_genotypes&, bool propagate);

  void parent_set_move(const individual_type& index, const individual_type& mate);
  void child_set_move (const individual_type& index, const individual_type& mate);
  void move_family    (const family_type& index);

  // Testing data
  bool informative_family(imodel&, const family_type& i) const;

  const subpedigree_type*  subpedigree;
  individual_const_pointer mother_ptr;
  individual_const_pointer father_ptr;
  child_vector             my_child_info;
  inconsistency_handler    err;
  fam_set                  set1, set2;

private:

   // Construction of calculation data structures.

   size_t build_family_data(imodel&, const family_type&);

   void build_child_data (imodel&, const individual_type& child_index, child_info_type&);

   // Removal of genotypes
   
   void remove_genotypes  (imodel&, const individual_type& par,
                           par_iterator begin, par_iterator end);

   bool remove_genotypes(imodel&, const child_info_type&);
};

// NOTE:  The bit class is used to store which genotypes of a particular
// individual are valid at any given time during the genotype elimination.
// This class is the bit_field.  While a vector<bool> may be substituted,
// several functions used (empty in particular) do not do the same things
// with these two classes.  Modification of the genotype elimination
// routines will be required.

#include "gelim/geno_eliminate.ipp"

}

#endif
