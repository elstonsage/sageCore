#ifndef MLM_CORRELATION_VERIFIER_H
#define MLM_CORRELATION_VERIFIER_H

#include "segreg/PedigreeDataSet.h"
#include "segreg/member_calculator.h"
#include "segreg/resid_sub_model.h"
#include "boost/array.hpp"

namespace SAGE {
namespace SEGREG {

/// \brief Checks the mlm correlations for bounds
///
/// The MlmCorrelationVerifier checks the constraints upon the mlm correlations
/// (actually adjustments).  These constraints are checked for each
/// correlated pair in the dataset.  Pair types that are checked are fs, ms, ss,
/// and fm.  See SEGREG documentation notes.
class MlmCorrelationVerifier
{
  public:
  
     MlmCorrelationVerifier(const PedigreeDataSet::SubpedigreeCursor& speds,
                            const residual_correlation_sub_model&     resid,
                            const binary_member_calculator&           bmc);

     bool current_correlation_estimates_are_valid() const;

     // Testing function not to be used in actual program
     void dump_valid_pairs() const;

  protected:

    typedef residual_correlation_sub_model::corr pair_class;
    typedef pair<FPED::MemberConstPointer,
                 FPED::MemberConstPointer>       pair_type;
    typedef list<pair_type>                      pair_list;
    
    bool is_valid_pair(const FPED::Member&, const FPED::Member&) const;

    void add_pair(const FPED::Member&, const FPED::Member&, pair_class);

    boost::array<pair_list, NUM_OF_CORRS> my_pairs_by_type;

    const binary_member_calculator&       my_bmc;
    const residual_correlation_sub_model& my_resids;
};

inline bool
  MlmCorrelationVerifier::is_valid_pair
    (const FPED::Member& i1, const FPED::Member& i2) const
{
  return my_bmc.is_member_valid(i1) && 
         my_bmc.is_member_valid(i2);
}

inline void
  MlmCorrelationVerifier::add_pair
    (const FPED::Member& i1, const FPED::Member& i2, pair_class p)
{
  my_pairs_by_type[p].push_back(make_pair(&i1, &i2));
}

}}

#endif

