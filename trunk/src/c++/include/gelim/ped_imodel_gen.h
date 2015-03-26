#ifndef PEDIMODEL_H
#define PEDIMODEL_H

//
// Pedigree Specific Inheritance Model Generation
//
// History: 0.1 gcw Initial Implementation
//
//
// Copyright (c) 2001 R. C. Elston.

#include <vector>
#include "mlocus/imodel.h"
#include "fped/fped.h"
#include "gelim/geno_eliminate.h"
#include "containers/bitfield.h"

namespace SAGE
{

/// \brief The pedigree_imodel_generator creates a new inheritance_model for the subpedigree it is given.
///
/// Given an inheritance_model and a subpedigree with phenotypes at that
/// model, it is often useful to be able to create an inheritance model
/// specific to that subpedigree.  This is because we can then isolate
/// those parts of the inheritance model that affect the subpedigree, and remove
/// everything else.  This includes alleles and phenotypes that are not
/// present, but more importantly, allows for the removal of penetrance
/// values that are invalid due to interactions between the subpedigree members
/// (genotype reduction and/or elimination).
///
/// The pedigree_imodel_generator creates a new inheritance_model for the
/// subpedigree it is given.  In the new inheritance_model, there is a 1
/// to 1 mapping between individual and their phenotype.  This new phentoype
/// has index one higher than the individual's in the subpedigree to
/// accomodate the missing phenotype.
///
/// CONSTRUCTION OPTIONS:
///
/// There are several boolean controls to specify how the imodel is to be
/// constructed.
/// 
/// - Remapping:
///
///   Remapping is the reduction of alleles that are not present into a single
///   'remap' allele. This has the advantage of making many algorithms faster
///   as only those alleles which might be present need be checked.  In general
///   these options should always be left as true in application code, but for
///   testing is allowed as an option.
/// 
///   There are two times remapping may be done.  The first is before genotype
///   elimination (if any) is performed.  This speeds up genotype elimination
///   at highly polymorphic loci.  The second is after genotype elimination is
///   performed and is for the removal of non-possible states that have
///   resulted because of genotype elimination.  In general, both should be
///   performed.
///
/// - Genotype Elimination:
///
///   Genotype elimination is the process of removing invalid genotypes from
///   individuals based upon the phenotypes of surrounding individuals. It is
///   especially useful for subpedigrees that contain many missing phenotype
///   individuals as often they can be restricted to a few valud states based
///   upon their relations.  More information can be found in the genotype
///   elimination header files.

class pedigree_imodel_generator
{
public:

  typedef FPED::Subpedigree                  subpedigree;
  typedef genotype_eliminator::handler       handler;

  pedigree_imodel_generator();
  
  /// @name Construction options
  //@{
  bool do_prior_remap() const;               ///< Do a remapping step prior to genotype elimination
  bool do_post_remap()  const;               ///< Do a remapping step prior to genotype elimination

  bool do_genotype_elimination() const;      ///< Do genotype elimination
  
  void set_prior_remap          (bool=true); ///< Set remapping step prior to G.E. to T/F
  void set_post_remap           (bool=true); ///< Set remapping step after to G.E. to T/F
  void set_genotype_elimination (bool=true); ///< Set genotype elimination to T/F
  //@}

  /// @name Construction/Processing
  //@{

  /// Given a subpedigree, create a subpedigree specific MLOCUS::inheritance_model,
  /// using remapping and genotype elimination if requested.  The base 
  /// MLOCUS::inheritance_model is taken from the multipedigree info based upon
  /// index \c m.
  ///
  /// The phenotype_id in this new model is the index of the individual
  /// in p + 1.  Phenotypes have internal names that should not be reported 
  /// outside of SAGE.  Returns an MLOCUS::inheritance_model if everything is fine,
  /// or an empty model if not (can be tested with genotype_informative() or 
  /// penetrance_informative())
  ///
  /// \param p The subpedigree in question
  /// \param m The index of the marker in the multipedigree info class
  ///
  /// \returns A new inheritance model based upon p and imod.
  MLOCUS::inheritance_model operator() (const subpedigree& p, size_t m) const;

  /// Given a subpedigree, create a subpedigree specific MLOCUS::inheritance_model,
  /// using remapping and genotype elimination if requested.  The inheritance model
  /// is taken from the function arguments, though the base MLOCUS::inheritance_model
  /// referred to by \c m may be queried for initial missing status and such.
  ///
  /// The phenotype_id in this new model is the index of the individual
  /// in p + 1.  Phenotypes have internal names that should not be reported 
  /// outside of SAGE.  Returns an MLOCUS::inheritance_model if everything is fine,
  /// or an empty model if not (can be tested with genotype_informative() or 
  /// penetrance_informative())
  ///
  /// \param p    The subpedigree in question
  /// \param m    The index of the marker in the multipedigree info class, used
  ///             for querying initial status information
  /// \param imod The inheritance model to use as a basis for the new inheritance
  ///             model.  Individual phenotypes of this model will be copied into
  ///             the new model.
  /// \param pids The phenotype ids in \c imod of the individuals in \c p
  ///
  /// \returns A new inheritance model based upon p and imod.
  MLOCUS::inheritance_model operator() (const subpedigree& p, size_t m,
                                        const MLOCUS::inheritance_model& imod,
                                        const vector<uint>& pids) const;
  //@}
  
  /** @name Additional inheritance_model information
   *    These functions provide additional information
   *    about the inheritance_model just returned. */
  //@{
  bool inconsistent() const;  ///< Were any inconsistencies found?
  bool informative()  const;  ///< Was the marker informative?

  //@}

  /// @name Output of errors
  //@{
  
  const handler& get_errors() const { return gelim.get_errors(); }

  void clear_errors() { gelim.clear_errors(); }

  //@}

protected:

  /// For each individual in the subpedigree which has a missing phenotype
  /// according to the pids vector, copy the missing phenotype from the original
  /// model into the new model.
  ///
  /// \param model The model to which we're copying.  It uses the
  ///              member's subpedigree index + 1 phenotype id schema.
  /// \param orig  The original model
  /// \param sped  The subpedigree
  /// \param pids  The member's phenotype ids in the original model.
  void copy_missing_members (MLOCUS::inheritance_model&       model,
                             const MLOCUS::inheritance_model& orig,
                             const subpedigree&               sped,
                             const vector<uint>&              pids) const;

  MLOCUS::inheritance_model empty_model(MLOCUS::inheritance_model, const subpedigree& sped)       const;
  MLOCUS::inheritance_model empty_model(const subpedigree&) const;

  /// Returns number of alleles remapped
  size_t do_remap(MLOCUS::inheritance_model&) const;

  /// Returns genotype elimination return value (see geno_eliminate.h)
  size_t do_gelim(const subpedigree&, MLOCUS::inheritance_model&, size_t marker) const;

  typedef std::map<uint, bool> RemapSetType;
  /// Find the set of alleles to be remapped
  void determine_allele_set(const MLOCUS::inheritance_model& i, RemapSetType&) const;

  bool prior_remap;               ///< Controls Prior Remapping
  bool post_remap;                ///< Controls Post Remapping
  bool geno_elim;                 ///< Controls Genotype Elimination

  mutable bool last_model_incon;  ///< Additional Information - last model was inconsistent
  mutable bool last_model_inform; ///< Additional Information - last model was informative

  mutable genotype_eliminator gelim;

};

#include "gelim/ped_imodel_gen.ipp"

}

#endif
