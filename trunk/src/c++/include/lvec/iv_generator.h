#ifndef __IV_GENERATOR
#define __IV_GENERATOR

//
//  The Inheritance Vector Generator - Generates inheritance vectors and
//      probabilities associated given a pedigree.
//
//  Author: Geoff Wedig
//
//  History:  0.1 Initial Implementation      Mar 4, 1998
//            1.0 Porting and Update          Apr    1998
//
//  Copyright (c) 1998 R. C. Elston
//

#include "boost/smart_ptr.hpp"
#include "gelim/ped_imodel_gen.h"
#include "lvec/inheritance_vector.h"

#define HFLAGS

namespace SAGE
{

class parent_vector;

/// a base class where iv_generator output is sent.
/** The iv_acceptor is a base class.  It accepts an inheritance_vector (or
 *  equivalence_class) and a probability associated with it and does some
 *  processing on it.  This processing may be generating a likelihood vector
 *  or pairs, or some other method.  
 * 
 *  The single function of accepting an equivalence class needs overloading.
 */
class iv_acceptor
{
public:

  // Required to make compiler happy.
  virtual ~iv_acceptor() { }
  
  typedef inheritance_vector::equivalence_class equivalence_class;

  void accept(equivalence_class e, double p) { accept(e, 0, p); }

  virtual void accept(equivalence_class e, equivalence_class dont_care_bits, double p) = 0;
};

/// Generates valid inheritance_vectors
/** The iv_generator takes an iv_acceptor at construction, Given a
 *  meiosis_map and an inheritance_model for a pedigree (must be the same)
 *  it generates the inheritance_vectors that are valid for this pedigree.
 *
 *  The iv_acceptor is passed and stored as a pointer.  The iv_generator should
 *  therefore never be used after it's iv_acceptor has gone out of scope or been
 *  deleted.  A new iv_acceptor should be given it instead.
 */
class iv_generator
{
public:

  typedef meiosis_map::subpedigree_const_pointer subpedigree_const_pointer;
  typedef meiosis_map::member_const_pointer      member_const_pointer;

  iv_generator() : iva() { }
  iv_generator(const boost::shared_ptr<iv_acceptor>& i) : iva() { set_iva(i); }
  
  ~iv_generator() { }
  
  /// @name Modification of constructor variables.
  //@{
  iv_acceptor&         get_iva() const { return *iva; }
  iv_acceptor&         set_iva(const boost::shared_ptr<iv_acceptor>& i)
  { iva = i; return *iva; }
  //@}
  
  /** In the build step, we assume that the inheritance_model passed is for
   * the pedigree section of the meiosis map and has had genotype
   * elimination and remapping done to it.  If the inheritance_model is
   * incorrect then the build will return false.  If elimination and
   * remapping has not been done, the algorithm may take much longer to
   * complete (but will complete, eventually) */
  bool build(const SAGE::meiosis_map&, const SAGE::MLOCUS::inheritance_model&,
             bool use_pop_freq = true) const;
  
protected:

  void build_ind(size_t, double, const SAGE::MLOCUS::inheritance_model&, 
                 parent_vector& pv) const;

  void evaluate(double prob) const;

#ifdef HFLAGS
  void set_bit(list<size_t>::iterator&, double prob)  const;
#endif

  boost::shared_ptr<iv_acceptor> iva;

  /// Processing Variables
  /** These variables are for internal processing, so they don't need to
   *  remain const.
   */
  //@{
  mutable SAGE::meiosis_map       mm;
  mutable inheritance_vector      iv;
  mutable vector<SAGE::MLOCUS::phased_genotype> genotypes;

#ifdef HFLAGS
  mutable list<size_t> hflags;  // homozygous flags
#endif
  //@}

  mutable bool use_pf;
};

}

#endif
