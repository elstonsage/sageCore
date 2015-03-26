//======================================================================
//
//  File:  Kernel.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//
//  See Kernel.h for summary information.
//
//======================================================================

#include "ageon/Kernel.h"

namespace SAGE {
namespace AO {

//======================================================================
//
//  Kernel(...)
//
//======================================================================
Kernel::Kernel(
  const Model      & mod, 
        const SAMPLING::PartitionedMemberDataSample & sample,
        int             t) 
  :
  my_analysis_type (t),
  my_mcc           (mod, sample, t),
  my_model         (mod), 
  my_sample        (sample)
{ }

//======================================================================
//
//  update()
//
//======================================================================
int
Kernel::update()
{
  int e = my_mcc.update();

  return e;
}

//======================================================================
//
//  get_mcc()
//
//======================================================================
const MemberCovariateCalculator &
Kernel::get_mcc() const
{
  return my_mcc;
}

//======================================================================
//
// get_sample_likelihood()
//
//======================================================================
double
Kernel::get_sample_likelihood()
{
  // 0. Set up local variables:

  log_double sample_likelihood (1.0),
             indiv_likelihood  (1.0);

  // 1. Loop through samplegroups and add up the likelihoods:

  for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
  {
    indiv_likelihood = get_indiv_likelihood(i);

    sample_likelihood *= indiv_likelihood;
  }

  // 3. Return likelihood:

  return sample_likelihood.get_log();
}

//======================================================================
//
//  get_indiv_likelihood(...)
//
//======================================================================
log_double
Kernel::get_indiv_likelihood(size_t indiv)
{
  // 0. Verify that the individual is valid before trying to calculate likelihood:

  if(!my_sample.isValid(indiv) || !my_sample.getMultipedigree().member_index(indiv).parent1())
    return log_double(1.0);

  // 1. There are three kinds of individual likelihoods, depending on the boolean 
  //    values AO_known, and affetedness.

  else
  {
    bool AO_known = my_sample.getField(my_model.AO_id).isAdjValuePresent(indiv),
         affected = my_sample.getAdjValue(indiv, my_model.aff_id);

         if(affected &&  AO_known) { return get_indiv_likelihood_affected_AO_known   (indiv); }
    else if(affected && !AO_known) { return get_indiv_likelihood_affected_AO_unknown (indiv); }
    else                           { return get_indiv_likelihood_unaffected          (indiv); }
  }
}

//======================================================================
//
//  get_indiv_likelihood_affected_AO_known(...)
//
//======================================================================
log_double 
Kernel::get_indiv_likelihood_affected_AO_known(size_t indiv)
{
  // 0. Set up local variables:

  size_t AO          = (size_t)my_sample.getAdjValue(indiv, my_model.AO_id);
  double suscept     = my_mcc   . get_genetic_suscept (indiv),
         mean        = my_mcc   . get_AO_mean         (indiv),
         stdev       = my_mcc   . get_AO_stdev        (indiv),
         AO_transf   = my_mcc   . get_AO_transf       (indiv),
         lambda1     = my_model . GetParameterMgr()(my_model.lambda1_mxid),
         lambda2     = my_model . GetParameterMgr()(my_model.lambda2_mxid),
         AO_term     = pow(AO + lambda2, lambda1 - 1.0) / stdev,
         pdf_term    = std_nor_density(AO_transf - mean, stdev);

  // 1. Calculate individual likelihood (allowing for truncation):

  log_double ind_lh(1.0);

  if(UseTruncation(my_analysis_type))
  {
    double zero_transf = my_mcc.get_zero_transf (),
           cdf_term    = normal_cdf(sign(lambda1) * (mean - zero_transf) / stdev);

    ind_lh = (suscept * pdf_term * AO_term) / cdf_term;
  }
  else
  {
    ind_lh = suscept * pdf_term * AO_term;
  }

  // 2. Return individual likelihood:

  return ind_lh;
}

//======================================================================
//
//  get_indiv_likelihood_affected_AO_unknown(...)
//
//======================================================================
log_double
Kernel::get_indiv_likelihood_affected_AO_unknown(size_t indiv)
{ 
  // 0. Set up local variables:
  
  double suscept     = my_mcc   . get_genetic_suscept (indiv),
         mean        = my_mcc   . get_AO_mean         (indiv),
         stdev       = my_mcc   . get_AO_stdev        (indiv),
         AE_transf   = my_mcc   . get_AE_transf       (indiv),
         cdf1        = normal_cdf((AE_transf - mean) / stdev);

  // 1. Calculate individual likelihood (allowing for truncation):

  log_double ind_lh(1.0);

  if(UseTruncation(my_analysis_type))
  {
    double zero_transf = my_mcc   . get_zero_transf (),
           lambda1     = my_model . GetParameterMgr ()(my_model.lambda1_mxid),
           cdf2        = normal_cdf(                (zero_transf - mean)        / stdev),
           cdf3        = normal_cdf(sign(lambda1) * (mean        - zero_transf) / stdev);

    ind_lh = (suscept * sign(lambda1) * (cdf1 - cdf2)) / cdf3;
  }
  else
  {
    ind_lh = suscept * cdf1;
  }

  // 2. Return individual likelihood:

  return ind_lh;
}

//======================================================================
//
//  get_indiv_likelihood_unaffected(...)
//
//======================================================================
log_double 
Kernel::get_indiv_likelihood_unaffected(size_t indiv)
{
  // 0. Set up local variables:
 
  log_double ind_lh (1.0),
             term1  (1.0),
             term2  (get_indiv_likelihood_affected_AO_unknown(indiv));

  // 1. Calculate likelihood:

  ind_lh = term1 - term2;
  
  // 2. Return likelihood:

  return ind_lh;
}

} // End namespace AO
} // End namespace SAGE
