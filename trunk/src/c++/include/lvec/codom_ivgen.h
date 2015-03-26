#ifndef __CODOM_IV_GENERATOR
#define __CODOM_IV_GENERATOR

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
#include "mlocus/imodel.h"
#include "fped/fped.h"
#include "lvec/inheritance_vector.h"
#include "lvec/iv_generator.h"

#undef HFLAGS

namespace SAGE
{

/** The codominant_iv_generator takes an iv_acceptor at construction, Given a
 *  meiosis_map and an inheritance model for a pedigree (must be the same or
 *  an error is called) it generates the inheritance vectors that are valid
 *  for this pedigree. It uses genotype elimination and marker remapping to
 *  be more efficient on each individual pedigree.
 */
class codominant_iv_generator
{
public:

  typedef FPED::SubpedigreeConstPointer ped_id;

  codominant_iv_generator() : my_iva() { }
  codominant_iv_generator(const boost::shared_ptr<iv_acceptor>& i) : my_iva() { set_iva(i); }
  
  ~codominant_iv_generator() { }
  
// Modification of constructor variables.

  iv_acceptor&         get_iva() const { return *my_iva; }
  iv_acceptor&         set_iva(const boost::shared_ptr<iv_acceptor>& i)
  { my_iva = i; return *my_iva; }

// And the build

  // Note.  We assume that the inheritance_model passed is codominant for
  // all typed individuals.  If the inheritance_model is not codominant,
  // the build will fail.
  bool build(const SAGE::meiosis_map&, const SAGE::MLOCUS::inheritance_model&,
             bool use_pop_freq = true) const;
  
protected:

  struct ind_type
  {
    long edge_number;
    long node1, node2;
    bool process;

    ind_type() : edge_number(-1), process(false) { }
  };

  void back(int) const;

  void build_ind(size_t i) const;
  void add_edge(int, int, int) const;

  void eval_ind(size_t, double) const;
  void evaluate(double) const;
  void set_bit(list<size_t>::const_iterator&, double)  const;

  boost::shared_ptr<iv_acceptor> my_iva;

// These variables are for internal processing, so they don't need to remain
// const.

  mutable vector<ind_type>             my_inds;
  mutable meiosis_map                  my_mmap;
  mutable SAGE::MLOCUS::inheritance_model            my_imodel;
  mutable inheritance_vector           my_ivector;
  mutable list<size_t>                    my_hflags;

  mutable double                       my_penetrance;

  mutable bool use_pf;
};

}

#endif
