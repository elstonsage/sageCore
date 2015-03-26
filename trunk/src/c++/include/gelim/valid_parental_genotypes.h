#ifndef VALID_PARENTAL_GENOTYPES_H
#define VALID_PARENTAL_GENOTYPES_H

#include <vector>
#include <utility>
#include "LSF/LSF.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "fped/fped.h"
#include "containers/bitfield.h"
#include "mlocus/imodel.h"
#include "gelim/inconsistency_handler.h"

#ifdef __KCC
using std::pair;
using std::vector;
#endif

namespace SAGE
{
/// Generates the valid genotypes and combinations for a pair of parents
/** Given a nuclear family, generate the set of valid parental genotypes
 *  based upon the children.  Additionally, generates a vector of valid and
 *  invalid genotypes for each of the parents (useful for genotype
 *  elimination and other algorithms)
 *
 *  This object is intended for use within the genotype_eliminator, as well
 *  as in other algorithms where a set of valid parental genotypes needs
 *  to be obtained.
 */
class valid_parental_genotypes
{
public:

  typedef FPED::Family                                  family_type;
  typedef MLOCUS::inheritance_model                     imodel;
  typedef pair<uint, uint>                              parental_genotype_pair;
  typedef vector<parental_genotype_pair>                parental_genotype_pair_vector;
  typedef vector<uint>                                  parental_genotype_vector;
  typedef parental_genotype_pair_vector::const_iterator genotype_pair_iterator;
  typedef parental_genotype_vector::const_iterator      genotype_iterator;

  valid_parental_genotypes();
  valid_parental_genotypes(const family_type&, const imodel&,
                           bool parents_phased = false);
  valid_parental_genotypes(const valid_parental_genotypes&);

  valid_parental_genotypes& operator=(const valid_parental_genotypes&);
  
  /// Calculation of valid pairs
  
  bool generate_valid_parental_genotypes(const family_type& fm,
                                         const imodel&      model,
                                         bool               parents_phased = false);

  /// Iteration over/Indexing the genotype pairs
  //@{  
  genotype_pair_iterator genotype_pair_begin() const;
  genotype_pair_iterator genotype_pair_end()   const;
    
  parental_genotype_pair genotype_pair(size_t i) const;

  size_t genotype_pair_count() const;
  //@}

  /// Iteration over the valid genotypes for parents
  //@{  
  genotype_iterator mother_valid_genotype_begin() const;
  genotype_iterator mother_valid_genotype_end()   const;
    
  uint mother_valid_genotype(size_t i) const;

  size_t mother_valid_genotype_count() const;

  genotype_iterator father_valid_genotype_begin() const;
  genotype_iterator father_valid_genotype_end()   const;
    
  uint father_valid_genotype(size_t i) const;

  size_t father_valid_genotype_count() const;
  //@}

  /// Iteration over the invalid genotypes for parents
  //@{  
  genotype_iterator mother_invalid_genotype_begin() const;
  genotype_iterator mother_invalid_genotype_end()   const;
    
  uint mother_invalid_genotype(size_t i) const;

  size_t mother_invalid_genotype_count() const;

  genotype_iterator father_invalid_genotype_begin() const;
  genotype_iterator father_invalid_genotype_end()   const;
    
  size_t father_invalid_genotype_count() const;

  uint father_invalid_genotype(size_t i) const;
  //@}

protected:

  void clean_internals();

  bool generate_valid_unphased_parental_genotypes
           (const family_type& fam,const imodel& model);
  bool generate_valid_phased_parental_genotypes
           (const family_type& fam,const imodel& model);

  parental_genotype_pair_vector my_genotype_pairs;

  parental_genotype_vector      my_mother_valid;
  parental_genotype_vector      my_father_valid;
  parental_genotype_vector      my_mother_invalid;
  parental_genotype_vector      my_father_invalid;
};

#include "gelim/valid_parental_genotypes.ipp"

}

#endif
