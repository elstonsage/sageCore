#include "freq/SingletonCalc.h"

namespace SAGE {
namespace FREQ {

//================================================================
//
//  Constructor
//
//================================================================
SingletonCalc::SingletonCalc(
  const MAXFUN::ParameterMgr   & mgr,
  const Sample                 & sample,
        bool                     use_inbreeding)
  : 
  my_sample         (sample),
  my_mgr            (mgr),
  my_geno_prob_calc (sample.getMarkerInfo().gmodel())
{ 
  my_geno_prob_calc.setUseInbreeding(use_inbreeding);
}

//================================================================
//
//  COPY Constructor
//
//================================================================
SingletonCalc::SingletonCalc(const SingletonCalc & other) :
  my_sample         (other.my_sample),
  my_mgr            (other.my_mgr),
  my_geno_prob_calc (other.my_geno_prob_calc)
{ }

//================================================================
//
//  calculateLikelihood()
//
//================================================================
log_double
SingletonCalc::calculateLikelihood()
{
  // Repopulate the genotype cache:
  my_geno_prob_calc.populateCache(my_mgr);

  // Calculate the likelihood:
  log_double lh(1.0);

  Sample::MemberVector::const_iterator member_end_itr = my_sample.getUnconnecteds().end();

  // Loop across all individuals:
  for(Sample::MemberVector::const_iterator ind = my_sample.getUnconnecteds().begin(); ind != member_end_itr; ++ind)
  {
    log_double ind_lh(0.0);

    const GenotypeInfoVector & genotypes = my_sample.getGenotypeInfoVector(**ind);
  
    GenotypeInfoVector::const_iterator genotype_end_itr = genotypes.end();

    for(GenotypeInfoVector::const_iterator g = genotypes.begin(); g != genotype_end_itr; ++g)
    {
      ind_lh += my_geno_prob_calc.getProb(g->genotype.get_id()) * g->penetrance;
    }

    lh *= ind_lh;
  }

  return lh;
}


} // End namespace FREQ
} // End namespace SAGE

