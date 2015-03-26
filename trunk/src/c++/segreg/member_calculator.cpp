//=========================================================================
//
//  File:	member_calculator
//
//  Purpose:	Calculates the crazy stuff.
//
//  Author:	Stephen Gross
//
//  History:	0.1  sag  Initial implementation	Jul 11 01
//              0.2  sag  Added documentation		Aug 13 01
//
//  Copyright (c) 2001 R. C. Elston
//=========================================================================

/// \file
///
/// The member calculators take care of calculating various values for each
/// member of each pedigree.  These values vary depending on the type of
/// model being fit (continuous, binary or onset) and sometimes on the
/// analysis type (FPMM, class A vs. class D, etc)

#include "segreg/member_calculator.h"

// This define is used for simplicity, and is undef'd at the end of the file
#define ASM ascertainment_sub_model

using namespace std;

bool finite2(double d) { return finite(d); }

namespace SAGE {
namespace SEGREG {

//lint --e{732}

//======================================================================
//                                                                     =
//  member_calculator_base() constructor                          =
//                                                                     =
//======================================================================

member_calculator_base::member_calculator_base(
  const FPED::Multipedigree & ped_data,
  const model               & mdl,
  bool                        use_asc)
  : my_ped_data (ped_data),
    mod         (mdl)
{
  // Deal with ascertainment

  if(use_asc && mod.ascer_sub_model.s_option() != ASM::none)
    use_ascertainment = true;
  else
    use_ascertainment = false;

  // Create the common set of data structures

  size_t total_inds   = 0;

  for(FPED::PedigreeConstIterator 
        pedigree_loop  = my_ped_data.pedigree_begin();
        pedigree_loop != my_ped_data.pedigree_end();
      ++pedigree_loop)
  {
    // Add a new member to pedigree_index_map to point to the data location

    //lint -e{534}
    pedigree_index_map.insert(make_pair(&*pedigree_loop, total_inds));

    total_inds += pedigree_loop->info().member_count();
  }

  last_pedigree_accessed = NULL;
  last_pedigree_location = (size_t) -1;

  member_classes.resize(total_inds, missing);
}

void
member_calculator_base::import_covariate_data
    (pedigree_const_pointer     ped,
     size_t                     member_index,
     size_t                     abs_mem_ref,
     vector<vector<double> >&   cov_data,
     const CovariateSubmodel&   cov_model)
{
  // Import the covariate data, terminating early if the
  // individual is invalidated.

  for(size_t i = 0; is_member_valid(abs_mem_ref) && i < cov_data.size(); ++i)
  {
    size_t trait_index = cov_model.covariates()[i].trait_index;

    //lint -e{534}
    import_trait(trait_index, ped, member_index, abs_mem_ref,
                 cov_data[i][abs_mem_ref]);
  }
}

//======================================================================
//                                                                     =
//  continuous_member_calculator() constructor                          =
//                                                                     =
//======================================================================

continuous_member_calculator::continuous_member_calculator(
  const FPED::Multipedigree & ped_data,
  const model               & mdl,
  bool                        use_asc)
  : member_calculator_base(ped_data, mdl, use_asc)
{
  // Check to see if the model is a continuous one.  If not, return without
  // doing anything

  if(mod.get_primary_trait_type() != pt_CONTINUOUS) return;

  // Import the data from the multipedigree and center it

  allocate_memory ();
  import_data     ();
  center_data     ();
}

//======================================================================
//                                                                     
//  calculate_composite_traits()
//                                                                     
//======================================================================

int
continuous_member_calculator::calculate_composite_traits()
{
  // Copy the primary traits into the composite_traits
  
  composite_traits = primary_traits; 

  // We only have to do this loop if there are covariates, so we check that first
  if(mod.comp_trait_sub_model.covariates().size())
  {
    for(size_t i = 0; i < member_classes.size(); ++i)
    {
      if(!is_member_valid(i)) continue; 

      double& composite_trait = composite_traits[i];

      for(size_t j = 0; j < mod.comp_trait_sub_model.covariates().size(); ++j)
        composite_trait += mod.comp_trait_sub_model.covariates()[j].coefficient *
                           comp_trait_data[j][i];
    }
  }
  
  // We only calculate the geometric mean if we are not using ascertainment
  // This assures that we don't recalculate the geometric mean for the
  // transform on just the ascertained subset.

  if(!use_ascertainment)
  {
    if(!mod.transf_sub_model.calculate_geom_mean(composite_traits))
      return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  // Transform the trait.

  if(!mod.transf_sub_model.transform(composite_traits))
    return segreg_errors::MCC_FAILED_TRANSFORM;

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     
//  calculate_expected_means()
//                                                                     
//======================================================================

int
continuous_member_calculator::calculate_expected_means()
{
  //lint --e{732}

  for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
  {
    double base_mean = mod.mean_sub_model.parameter(genotype);
    
    // This is somewhat clever.  For each value in primary trait, if it is
    // finite (ie, not QNAN, and the member is valid), it puts the expected
    // mean in the expected means vector, and otherwise puts the value from
    // primary_traits (QNAN, since the pt cannot be infinite)

    //lint -e{534} <- Ignoring return
    replace_copy_if(primary_traits.begin(), primary_traits.end(),
                    expected_means[genotype].begin(), finite2, base_mean);

    // We only have to do this loop if there are covariates, so we check
    // that first

    if(mod.mean_cov_sub_model.covariates().size())
    {
      for(size_t i = 0; i < member_classes.size(); ++i)
      {
        if(!is_member_valid(i)) continue;

        double& expected_mean = expected_means[genotype][i];

        for(size_t j = 0; j < mod.mean_cov_sub_model.covariates().size(); ++j)
        {
          expected_mean += mod.mean_cov_sub_model.covariates()[j].coefficient *
                           mean_cov_data[j][i];

          if(mod.mean_cov_sub_model.covariates()[j].has_interaction)
          {
            double tau = mod.mean_cov_sub_model.covariates()[j].i_taus[genotype];

            expected_mean += mod.mean_sub_model.parameter(genotype) *
                             mean_cov_data[j][i] * tau;
          }
        }
      }
    }

    if(!mod.transf_sub_model.transform(expected_means[genotype]))
      return segreg_errors::MCC_FAILED_TRANSFORM;
  }
  
  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     
//  calculate_expected_polygenic_means()
//                                                                     
//======================================================================

int
continuous_member_calculator::calculate_expected_polygenic_means()
{
  for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
  {
    for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
    {
      double polygenic_adjustment = mod.fpmm_sub_model.mean(k);

      // Copy the expected_means into the polygenic_means
      expected_polygenic_means[k][genotype] = expected_means[genotype];
      
      // Adjust the values by the polygenic adjustment
      
      //lint -e{534} <- Ignoring return
      transform(expected_polygenic_means[k][genotype].begin(),
                expected_polygenic_means[k][genotype].end(),
                expected_polygenic_means[k][genotype].begin(),
                bind1st(plus<double>(), polygenic_adjustment));
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calc_exp_var(...)                                                 =
//                                                                     =
//======================================================================
void
continuous_member_calculator::calc_exp_var(
        size_t abs_mem_ref,
        genotype_index genotype)
{
  double& expected_var = expected_variances[genotype][abs_mem_ref];

  double  tau = 0.0;

  for(size_t i = 0; i < mod.var_cov_sub_model.covariates().size(); ++i)
  {
    expected_var += mod.var_cov_sub_model.covariates()[i].coefficient *
                    var_cov_data[i][abs_mem_ref];

    if(mod.var_cov_sub_model.covariates()[i].has_interaction)
    {
      tau = mod.var_cov_sub_model.covariates()[i].i_taus[genotype];

      expected_var += mod.var_sub_model.parameter(genotype) *
                      var_cov_data[i][abs_mem_ref] * tau;
    }
  }
}
//======================================================================
//                                                                     
//  calculate_expected_variances()
//                                                                     
//======================================================================

int
continuous_member_calculator::calculate_expected_variances()
{
  for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
  {
    double expected_var = mod.var_sub_model.parameter(genotype);

    // This is somewhat clever.  For each value in primary trait, if it is
    // finite (ie, not QNAN, and the member is valid), it puts the expected
    // mean in the expected means vector, and otherwise puts the value from
    // primary_traits (QNAN, since the pt cannot be infinite)

    //lint -e{534}
    replace_copy_if(primary_traits.begin(), primary_traits.end(),
                    expected_variances[genotype].begin(), finite2, expected_var);

    // We only have to do this loop if there are covariates, so we check
    // that first

    if(mod.var_cov_sub_model.covariates().size())
    {
      for(size_t i = 0; i < member_classes.size(); ++i)
      {
        if(is_member_valid(i))
        {
          calc_exp_var(i, genotype);

          if(get_expected_variance(i, genotype) <= 0.0)
            return segreg_errors::MCC_VARIANCE_INVALID;
        }
      }
    }
  }
  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calc_exp_sd(...)                                                   =
//                                                                     =
//======================================================================
void
continuous_member_calculator::calc_exp_sd(
	size_t abs_mem_ref,
	genotype_index genotype)
{
  expected_sds[genotype][abs_mem_ref] = 
    sqrt(get_expected_variance(abs_mem_ref,genotype));
}

//======================================================================
//                                                                     
//  calculate_expected_sds()
//                                                                     
//======================================================================

int
continuous_member_calculator::calculate_expected_sds()
{
  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(is_member_valid(i))
      for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
        calc_exp_sd(i,genotype);
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calculate_standardizations()                                       =
//                                                                     =
//======================================================================

int
continuous_member_calculator::calculate_standardizations()
{
  // Calculate the new standardizations.

  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    calc_stan(i,index_AA);
    calc_stan(i,index_AB);
    calc_stan(i,index_BB);
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calculate_ascertained_standardizations()                           =
//                                                                     =
//======================================================================

int
continuous_member_calculator::calculate_ascertained_standardizations()
{
  // Determine Thigh and Tlow, if needed

  double Thigh = mod.ascer_sub_model.thresh_high();
  double Tlow  = mod.ascer_sub_model.thresh_low();

  // For each non-nan value, we have to transform it.

  if(!SAGE::isnan(Thigh))
  {
    if(!mod.transf_sub_model.transform(Thigh))
      return segreg_errors::MCC_FAILED_TRANSFORM;
  }
  if(!SAGE::isnan(Tlow))
  {
    if(!mod.transf_sub_model.transform(Tlow))
      return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  // Calculate the new ascertained standardizations.

  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    calc_asc_stan(i,index_AA, Thigh, Tlow);
    calc_asc_stan(i,index_AB, Thigh, Tlow);
    calc_asc_stan(i,index_BB, Thigh, Tlow);
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calculate_polygenic_standardizations()                             =
//                                                                     =
//======================================================================

int
continuous_member_calculator::calculate_polygenic_standardizations()
{
  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
    {
      calc_pg_stan(i,index_AA,k);
      calc_pg_stan(i,index_AB,k);
      calc_pg_stan(i,index_BB,k);
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calculate_estimated_standardizations()                             =
//                                                                     =
//======================================================================

int
continuous_member_calculator::calculate_estimated_standardizations()
{
  // Calculate the new est. standardizations.

  vector<double>          tempv1(3);
  vector<vector<double> > tempv2(3, tempv1);

  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;
  
    for(genotype_index genotype_mother  = index_AA;
                       genotype_mother != index_INVALID;
                     ++genotype_mother)
    {
      for(genotype_index genotype_father  = index_AA;
                         genotype_father != index_INVALID;
                       ++genotype_father)
      {
        calc_est_stan(i, genotype_mother, genotype_father);

        // If we get a NaN estimated standardization, the variance is too
        // small, causing the estimated standardization to not be
        // calculatable.

        if(SAGE::isnan(get_estimated_standardization(i, genotype_mother, genotype_father)))
          return segreg_errors::MCC_VARIANCE_INVALID;
      }
    }
  }

  return segreg_errors::EVAL_OK;
}

void
continuous_member_calculator::allocate_memory()
{
  // Create some temporaries.

  vector<double>          temp2    (member_classes.size(), 0);
  vector<vector<double> > temp5    (3, temp2);    // Creates a 3 x i

  // Allocate primary storage

  primary_traits = temp2;

  comp_trait_data .resize(mod.comp_trait_sub_model.get_covariate_count(), temp2);
  mean_cov_data   .resize(mod.mean_cov_sub_model  .get_covariate_count(), temp2);
  var_cov_data    .resize(mod.var_cov_sub_model   .get_covariate_count(), temp2);

  // Allocate storage of calculated values.

  composite_traits             = temp2;
  expected_means               = temp5;
  expected_variances           = temp5;
  expected_sds                 = temp5;

  standardizations             = temp5;
  ascertained_standardizations = temp5;

  if(mod.get_model_class() == model_D)
  {
    estimated_standardizations .resize (3, temp5);
  }
  else if(mod.get_model_class() == model_FPMM)
  {
    expected_polygenic_means   .resize (mod.fpmm_sub_model.max_pgt(), temp5);
    polygenic_standardizations .resize (mod.fpmm_sub_model.max_pgt(), temp5);
  }
}

//======================================================================
//                                                                     =
//  import_data()                                                       =
//                                                                     =
//======================================================================

void
continuous_member_calculator::import_data()
{
  // Determine our trait index for the primary trait

  const string& pr_trait_name  = mod.get_primary_trait();
  size_t        pr_trait_index = my_ped_data.info().trait_find(pr_trait_name);

  size_t abs_mem_ref = 0;

  for(FPED::PedigreeConstIterator 
        pedigree_loop  = my_ped_data.pedigree_begin();
        pedigree_loop != my_ped_data.pedigree_end();
      ++pedigree_loop)
  {
    for(size_t member_loop = 0;
               member_loop < pedigree_loop->info().member_count();
             ++member_loop, ++abs_mem_ref)
    {
      // Assign new member.  Note that since we don't know the member's data
      // yet, we simply set to 'actual', assuming that the member will be
      // valid and no ascertainment weirdness is in effect.  This will be
      // corrected later.

      member_classes[abs_mem_ref] = actual;

      FPED::MemberConstPointer mem = &pedigree_loop->member_index(member_loop);

      // Only import all the data if we're not using ascertainment or the
      // individual is in C
      if(!use_ascertainment || mod.ascer_sub_model.is_ind_in_C(mem))
      {
        //lint -e{534}
        import_trait(pr_trait_index,   &*pedigree_loop, member_loop, abs_mem_ref,
                     primary_traits [abs_mem_ref]);

        // Import the covariate data, terminating early if the
        // individual is invalidated.

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              comp_trait_data, mod.comp_trait_sub_model);

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              mean_cov_data, mod.mean_cov_sub_model);


        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              var_cov_data, mod.var_cov_sub_model);

        if(use_ascertainment && is_member_valid(abs_mem_ref))
        {
          switch(mod.ascer_sub_model.get_ind_type(mem))
          {
            case ASM::not_specified: member_classes[abs_mem_ref] = missing;    break;
            case ASM::actual:        member_classes[abs_mem_ref] = actual;     break;
            case ASM::gte_thresh:    member_classes[abs_mem_ref] = gte_thresh; break;
            case ASM::lte_thresh:    member_classes[abs_mem_ref] = lte_thresh; break;

            case ASM::thresh_indic:
            case ASM::onset:         SAGE_internal_error(); // These should never occur!

          }
        }
      }
      else
        member_classes[abs_mem_ref] = missing;
        
      // Test the member for validity.  If invalid, set the calculated
      // values of the trait, mean, variance and standard deviation to QNAN. 
      // This saves time later, since we never have to calculate them.
      if(!is_member_valid(abs_mem_ref))
      {
        primary_traits [abs_mem_ref] = QNAN;

        for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
        {
          expected_means     [genotype][abs_mem_ref] = QNAN;
          expected_variances [genotype][abs_mem_ref] = QNAN;
          expected_sds       [genotype][abs_mem_ref] = QNAN;
          
          if(mod.get_model_class() == model_FPMM)
          {
            for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
              expected_polygenic_means[k][genotype][abs_mem_ref] = QNAN;
          }
        }
      }
    } 
  } 
}
 
//======================================================================
//                                                                     =
//  center_data()                                                      =
//                                                                     =
//======================================================================
void
continuous_member_calculator::center_data()
{
  // Calculate all the means -- This is *very* non-optimal, but I can't get
  // around it right now.  FIXME! XXX GCW 2002-03-28

  model* m = const_cast<model*>(&mod);

  m->comp_trait_sub_model.calculate_covariate_means(my_ped_data);
  m->mean_cov_sub_model  .calculate_covariate_means(my_ped_data);
  m->var_cov_sub_model   .calculate_covariate_means(my_ped_data);

  // Composite covariates:

  for(size_t trait_num = 0; trait_num < mod.comp_trait_sub_model.covariates().size(); ++trait_num)
  {
    double mean = mod.comp_trait_sub_model.get_covariate_mean(trait_num);

    center_trait(comp_trait_data[trait_num], mean);
  }

  // Mean covariates:

  for(size_t trait_num = 0; trait_num < mod.mean_cov_sub_model.covariates().size(); ++trait_num)
  {
    double mean = mod.mean_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(mean_cov_data[trait_num], mean);
  }

  // Variance covariates:

  for(size_t trait_num = 0; trait_num < mod.var_cov_sub_model.covariates().size(); ++trait_num)
  {
    double mean = mod.var_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(var_cov_data[trait_num], mean);
  }
}

//======================================================================
//                                                                     =
//  calc_stan(...)                                      =
//                                                                     =
//======================================================================

void
continuous_member_calculator::calc_stan(
    size_t abs_mem_ref,genotype_index genotype)
{
  //=========================================//
  //Standardization is given by the equation://
  //                                         //
  //        [ t(i) - O(u,i) ]                //
  //      -----------------------            //
  //              N(u,i)                     //
  //                                         //
  //=========================================//

  double& stan = standardizations[genotype][abs_mem_ref];

  double comp_trait      = get_composite_trait (abs_mem_ref          );
  double exp_mean        = get_expected_mean   (abs_mem_ref, genotype);
  double exp_sd          = get_expected_sd     (abs_mem_ref, genotype);

  stan = (comp_trait - exp_mean) / exp_sd;
}

//======================================================================
//                                                                     =
//  calc_asc_stan(...)                                                 =
//                                                                     =
//======================================================================
/// Ascertained Standardization is given by the equation:
///  
///  \f[
///     \frac{ t' - \Theta_u(i) } 
///          {      \eta_u(i)   }
///  \f]
///   
///  where:
///  
///  \f[
///     t' =  \left\{ \begin{array}{ll}
///             0        & \mbox{ if member is of type {\tt missing     } }
///         \\  t_i      & \mbox{ if member is of type {\tt actual      } }
///         \\  T_{high} & \mbox{ if member is of type {\tt gte\_thresh } }
///         \\  T_{low}  & \mbox{ if member is of type {\tt lte\_thresh } }
///                   \end{array}
///           \right.
///  \f]

// } <=== done because of the open bracket in formula above to not mix up
//        the editor when completing brackets.

void
continuous_member_calculator::calc_asc_stan
    (size_t abs_mem_ref,genotype_index genotype, double Th, double Tl)
{
  //Standardization is given by the equation:
  //
  //        [ t' - O(u,i) ] 
  //      -----------------------
  //              N(u,i)
  // 
  // where:
  //
  // t' =  0     if member is of type missing
  //    =  ti    if member is of type actual
  //    =  Thigh if member is of type gte_thresh
  //    =  Tlow  if member is of type lte_thresh
  //

  double& stan = ascertained_standardizations[genotype][abs_mem_ref];

  double comp_trait = 0.0;

  switch(member_classes[abs_mem_ref])
  {
    case missing :
      comp_trait = 0;
      break;
    
    case actual :
      comp_trait = get_composite_trait (abs_mem_ref);
      break;
    
    case gte_thresh :
      comp_trait = Th;
      break;

    case lte_thresh :
      comp_trait = Tl;
      break;

    case age_at_onset : // Never happens with continuous
      break;
  }

  double exp_mean = get_expected_mean   (abs_mem_ref, genotype);
  double exp_sd   = get_expected_sd     (abs_mem_ref, genotype);

  stan = (comp_trait - exp_mean) / exp_sd;
}

//======================================================================
//                                                                     =
//  calc_pg_stan(...)                                                  =
//                                                                     =
//======================================================================

void
continuous_member_calculator::calc_pg_stan(
        size_t         abs_mem_ref,
        genotype_index genotype,
        size_t         polygenotype)
{
  //=========================================//
  //Standardization is given by the equation://
  //                                         //
  //        [ t(i) - O(u,v,i) ]              //
  //      -----------------------            //
  //              N(u,i)                     // Please note that this equation is
  //                                         // under review regarding the use of
  //=========================================// genotype-specific standard deviation

  double& stan  = polygenic_standardizations[polygenotype][genotype][abs_mem_ref];

  double comp_trait  = get_composite_trait (abs_mem_ref                      );
  double exp_mean    = get_expected_mean   (abs_mem_ref,genotype,polygenotype);
  double exp_sd      = get_expected_sd     (abs_mem_ref,genotype             );

  stan = (comp_trait - exp_mean) / exp_sd;
}

//======================================================================
//                                                                     =
//  calc_est_stan(...)                                                 =
//                                                                     =
//======================================================================
void
continuous_member_calculator::calc_est_stan(
           size_t         abs_mem_ref,
           genotype_index genotype_mother,
           genotype_index genotype_father)
{
  double sigma_hat                 = 0.0,
         variance_hat              = 0.0,
         sd_hat                    = 0.0,
         upsilon                   = 0.0,
         prob                      = 0.0,
         power                     = 0.0,
         numerator_sigma           = 0.0,
         numerator_var             = 0.0,
         denominator               = 0.0;

  double& stan  = estimated_standardizations[genotype_mother][genotype_father][abs_mem_ref];

  for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
  {
    prob             = mod.transm_sub_model.prob(genotype,genotype_mother,genotype_father);
    power            = get_standardization(abs_mem_ref,genotype);
    upsilon          = prob * exp(-power*power/2.0);
    numerator_sigma += upsilon * get_expected_mean    (abs_mem_ref,genotype); 
    numerator_var   += upsilon * get_expected_variance(abs_mem_ref,genotype);
    denominator     += upsilon;
  }

  // If the denominator is 0, then the standardization cannot be calculated. 
  // This should only happen when the variance is small relative to the
  // mean, causing the standardization to be large, and the exp term (see
  // above) to be smaller than the smallest representable double.

  if(denominator == 0.0)
  {
    stan = QNAN;

    return;
  }

  sigma_hat    = numerator_sigma / denominator;
  variance_hat = numerator_var   / denominator;
  sd_hat       = sqrt(variance_hat);

  stan = (get_composite_trait(abs_mem_ref) - sigma_hat) / sd_hat;
}  

//======================================================================
//                                                                     =
//  binary_member_calculator() constructor                          =
//                                                                     =
//======================================================================

binary_member_calculator::binary_member_calculator(
  const FPED::Multipedigree & ped_data,
  const model               & mdl,
  bool                        use_asc)
  : member_calculator_base(ped_data, mdl, use_asc)
{
  // Check to see if the model is a binary one.  If not, return without
  // doing anything

  if(mod.get_primary_trait_type() != pt_BINARY) return;

  // Allocate memory for all our data, import the data and center it

  allocate_memory();
  import_data();
  center_data();
}

void binary_member_calculator::allocate_memory()
{
  // Resize all our vectors

  vector<double>          temp2    (member_classes.size(), QNAN);
  vector<vector<double> > temp5    (3, temp2); // genotype x inds

  // Allocate memory for each of our primary data variables.  We set
  // them to 'invalid' values for safety.

  affections    .resize(member_classes.size(),                        false);
  susc_cov_data .resize(mod.susc_cov_sub_model.get_covariate_count(), temp2);

  // Resize the calculation vectors

  expected_susceptibilities = temp5;
  penetrances               = temp5;

  if(mod.get_model_class() == model_FPMM)
    polygenic_penetrances   .resize (mod.fpmm_sub_model.max_pgt(), temp5);
}

//======================================================================
//                                                                     
//  calculate_expected_susceptibilities()
//                                                                     
//======================================================================

int
binary_member_calculator::calculate_susceptibilities()
{
  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
    {
      expected_susceptibilities[genotype][i] = calc_exp_susc(i, genotype);
    }
  }

  return segreg_errors::EVAL_OK;
}

int
binary_member_calculator::calculate_penetrances()
{
  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    double aff  = int_get_aff_status(i);

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
    {
      double mean = int_get_expected_susc(i, genotype);

      penetrances[genotype][i] = calculate_penetrance(mean, aff);
    }
  }

  return segreg_errors::EVAL_OK;
}

int
binary_member_calculator::calculate_polygenic_penetrances()
{
  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    double aff  = int_get_aff_status(i);

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
    {
      double mean = int_get_expected_susc(i, genotype);
      
      for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
      {
        double poly_mean = mean + mod.fpmm_sub_model.mean(k);

        polygenic_penetrances[k][genotype][i] = calculate_penetrance(poly_mean, aff);
      }
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  import_data()                                                       =
//                                                                     =
//======================================================================

void
binary_member_calculator::import_data()
{
  // Determine our trait index for affection

  const string& aff_trait_name  = mod.get_primary_trait();
  size_t        aff_trait_index = my_ped_data.info().trait_find(aff_trait_name);

  size_t abs_mem_ref = 0;

  for(FPED::PedigreeConstIterator 
        pedigree_loop  = my_ped_data.pedigree_begin();
        pedigree_loop != my_ped_data.pedigree_end();
      ++pedigree_loop)
  {
    // Add a new member to pedigree_index_map to point to the data location

    //lint -e{534}
    pedigree_index_map.insert(make_pair(&*pedigree_loop,abs_mem_ref));

    for(size_t member_loop = 0;
               member_loop < pedigree_loop->info().member_count();
             ++member_loop, ++abs_mem_ref)
    {
      // Assign a new member.  Note that since we don't know the member's
      // data yet, we simply set to available, assuming that the member will
      // be valid.

      member_classes[abs_mem_ref] = available;

      // Get the member pointer.

      FPED::MemberConstPointer mem = &pedigree_loop->member_index(member_loop);

      // Only import all the data if we're not using ascertainment or the
      // individual is in C
      if(!use_ascertainment || mod.ascer_sub_model.is_ind_in_C(mem))
      {
        //lint -e{534}
        import_trait(aff_trait_index,   &*pedigree_loop, member_loop, abs_mem_ref,
                     affections [abs_mem_ref]);

        // Import the susceptibility data, terminating early if the
        // individual is invalidated.

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              susc_cov_data, mod.susc_cov_sub_model);
      }
      else
        member_classes[abs_mem_ref] = missing;

      // Test the member for validity.  If invalid, set the penetrances to
      // 1.0.  This saves time later, since we never have to calculate them.
      if(!is_member_valid(abs_mem_ref))
      {
        penetrances[0][abs_mem_ref] = 1.0;
        penetrances[1][abs_mem_ref] = 1.0;
        penetrances[2][abs_mem_ref] = 1.0;

        if(mod.get_model_class() == model_FPMM)
        {
          for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
          {
            polygenic_penetrances[k][0][abs_mem_ref] = 1.0;
            polygenic_penetrances[k][1][abs_mem_ref] = 1.0;
            polygenic_penetrances[k][2][abs_mem_ref] = 1.0;
          }
        }
      }
    } 
  } 
}
 
//======================================================================
//                                                                     =
//  center_data()                                                      =
//                                                                     =
//======================================================================
void
binary_member_calculator::center_data()
{
  // Calculate all the means -- This is *very* non-optimal (the const cast),
  // but I can't get around it right now.  FIXME! XXX GCW 2002-03-28

  model* m = const_cast<model*>(&mod);

  m->susc_cov_sub_model  .calculate_covariate_means(my_ped_data);

  // Susceptibility covariates:

  for(size_t trait_num = 0; trait_num < mod.susc_cov_sub_model.covariates().size(); ++trait_num)
  {
    double mean = mod.susc_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(susc_cov_data[trait_num], mean);
  }
}

//======================================================================
//                                                                     =
//  calc_exp_susc(...)                                                 =
//                                                                     =
//======================================================================
double
binary_member_calculator::calc_exp_susc(
        size_t abs_mem_ref,
        genotype_index genotype)
{
  double expected_susceptibility = 0.0,
         tau                     = 0.0;

  if(!is_member_valid(abs_mem_ref))
  {
    expected_susceptibility = numeric_limits<double>::quiet_NaN();
  }
  else
  {
    expected_susceptibility += mod.susc_sub_model.parameter(genotype);

    for(size_t i = 0; i < mod.susc_cov_sub_model.covariates().size(); ++i)
    {
      expected_susceptibility += mod.susc_cov_sub_model.covariates()[i].coefficient *
                                 susc_cov_data[i][abs_mem_ref];

      if(mod.susc_cov_sub_model.covariates()[i].has_interaction)
      {
        tau = mod.susc_cov_sub_model.covariates()[i].i_taus[genotype];

        expected_susceptibility += mod.susc_sub_model.parameter(genotype) *
                                   susc_cov_data[i][abs_mem_ref] * tau;
      }
    }
  }

  return expected_susceptibility;
}

//=====================================================================================
//=====================================================================================
//
//  onset_member_calculator
//
//=====================================================================================
//=====================================================================================

onset_member_calculator::onset_member_calculator(
  const FPED::Multipedigree & ped_data,
  const model               & mdl,
  bool                        use_asc)
  : member_calculator_base(ped_data, mdl, use_asc)
{
  // Check to see if the model is a onset one.  If not, return without
  // doing anything

  if(mod.get_primary_trait_type() != pt_ONSET) return;

  // Allocate memory for all our data.

  allocate_memory();

  // Import the data from the multipedigree

  import_data();

  // Calculate our geometric mean.  If we can't do this, there's no point in
  // going on.  We only do this if we're not ascertained, as we want to use
  // the whole sample otherwise.  If we are ascertained, we assume that the
  // non-ascertained onset calculator has done its own job.

  if(!use_ascertainment)
  {
    mod.transf_sub_model.set_clear_each_sync(false);

    valid_geom_mean = mod.transf_sub_model.calculate_geom_mean(age_onsets);

    if(!valid_geom_mean)
      return;
  }
  else
    valid_geom_mean=true;

  // Center the data

  center_data();

}

void onset_member_calculator::allocate_memory()
{
  // Resize all our vectors

  vector<double>          temp2    (member_classes.size(), QNAN);
  vector<vector<double> > temp5    (3, temp2);

  // We need to allocate memory for each of our primary variables.  We set
  // them to 'invalid' values for safety.

  affections.resize(member_classes.size(), false);

  age_onsets = temp2;
  age_exams  = temp2;

  mean_cov_data .resize(mod.mean_cov_sub_model.get_covariate_count(), temp2);
  var_cov_data  .resize(mod.var_cov_sub_model .get_covariate_count(), temp2);
  susc_cov_data .resize(mod.susc_cov_sub_model.get_covariate_count(), temp2);

  // Resize the expected means correctly

  size_t exp_size = 0;

  if(mod.ons_sub_model.m_option() == onset_sub_model::m_A)
    exp_size = mod.fpmm_sub_model.max_pgt();
  else
    exp_size = 1;

  expected_means.resize  (exp_size, temp5);

  expected_alphas = temp5;

  // Resize the expected susceptibilities correctly

  if(mod.ons_sub_model.m_option() == onset_sub_model::m_S)
    exp_size = mod.fpmm_sub_model.max_pgt();
  else
    exp_size = 1;

  expected_susceptibilities.resize (exp_size, temp5);

  // Resize the temporary vector

  temp_calc_space = temp5;
}

//======================================================================
//                                                                     =
//  import_data()                                                       =
//                                                                     =
//======================================================================

void
onset_member_calculator::import_data()
{
  // Determine if we're ascertained and by_onset.  This has effects later
  // when classifying individuals

  bool by_onset = false;

  if(use_ascertainment &&
     mod.ascer_sub_model.v_option() == ASM::onset)
  {
    by_onset = true;
  }

  // Determine our trait indices for affection, age_onset and age_exam

  const string& aff_trait_name  = mod.ons_sub_model.affection_status();
  size_t        aff_trait_index = my_ped_data.info().trait_find(aff_trait_name);

  const string& onset_trait_name  = mod.ons_sub_model.age_of_onset();
  size_t        onset_trait_index = my_ped_data.info().trait_find(onset_trait_name);

  const string& exam_trait_name  = mod.ons_sub_model.age_at_exam();
  size_t        exam_trait_index = my_ped_data.info().trait_find(exam_trait_name);

  // Grab all the values we can.

  size_t abs_mem_ref = 0;

  for(FPED::PedigreeConstIterator 
        pedigree_loop  = my_ped_data.pedigree_begin();
        pedigree_loop != my_ped_data.pedigree_end();
      ++pedigree_loop)
  {
    // Add a new member to pedigree_index_map to point to the data location

    //lint -e{534}
    pedigree_index_map.insert(make_pair(&*pedigree_loop,abs_mem_ref));

    for(size_t member_loop = 0;
               member_loop < pedigree_loop->member_count();
             ++member_loop, ++abs_mem_ref)
    {
      // Assign a new member.  Note that since we don't know the member's
      // data yet, we simply set to age_of_onset, though that's likely to be
      // incorrect.

      member_classes[abs_mem_ref] = age_of_onset;

      // Get the member pointer.

      FPED::MemberConstPointer mem = &pedigree_loop->member_index(member_loop);

      // Only import all the data if we're not using ascertainment or the
      // individual is in C
      if(!use_ascertainment || mod.ascer_sub_model.is_ind_in_C(mem))
      {
        bool af = import_trait(aff_trait_index,   &*pedigree_loop, member_loop,
                               abs_mem_ref, affections [abs_mem_ref]);
        bool ao = import_trait(onset_trait_index, &*pedigree_loop, member_loop,
                               abs_mem_ref, age_onsets [abs_mem_ref]);
        bool ax = import_trait(exam_trait_index,  &*pedigree_loop, member_loop,
                               abs_mem_ref, age_exams  [abs_mem_ref]);

        // Determine individual classification.

        // Classification goes like this:
        //
        // If we have an affection status:
        //   If we have an age of onset AND affection is true:
        //     Clear age of exam (otherwise might be transformed)
        //     If !by onset -> age_of_onset
        //     else         -> age_at_onset
        //   else if we have age at exam:
        //     Clear age of onset (otherwise might be transformed)
        //     If affected  -> age_at_exam_aff
        //     else         -> age_at_exam_unaff
        //   else           -> invalid (missing)
        // else             -> invalid (missing)

        if(af)
        {
          if(ao && affections[abs_mem_ref])
          {
            age_exams[abs_mem_ref] = QNAN;

            if(!by_onset)
              member_classes[abs_mem_ref] = age_of_onset;
            else
              member_classes[abs_mem_ref] = age_at_onset;
          }
          else if(ax)
          {
            age_onsets[abs_mem_ref] = QNAN;

            if(affections[abs_mem_ref])
              member_classes[abs_mem_ref] = age_at_exam_aff;
            else
              member_classes[abs_mem_ref] = age_at_exam_unaff;
          }
          else
            member_classes[abs_mem_ref] = missing;
        }
        else
        {
          member_classes[abs_mem_ref] = missing;
        }

        // We continue importing traits only if the member continues to be valid

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              mean_cov_data, mod.mean_cov_sub_model);

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              var_cov_data, mod.var_cov_sub_model);

        import_covariate_data(&*pedigree_loop, member_loop, abs_mem_ref,
                              susc_cov_data, mod.susc_cov_sub_model);

        // If we've got an invalid member, we reset their age of onset and exam
        // so they don't get used for the transformation.

      }
      else
        member_classes[abs_mem_ref] = missing;

      if(!is_member_valid(abs_mem_ref))
      {
        age_onsets[abs_mem_ref] = QNAN;
        age_exams [abs_mem_ref] = QNAN;

        for(genotype_index geno = index_AA; geno != index_INVALID; ++geno)
        {
          expected_alphas[geno][abs_mem_ref] = QNAN;

          temp_calc_space[geno][abs_mem_ref] = QNAN;
        }
      } 
    }
  } 
}

//======================================================================
//                                                                     =
//  center_data()                                                      =
//                                                                     =
//======================================================================
void
onset_member_calculator::center_data()
{
  // Calculate all the means -- This is *very* non-optimal, but I can't get
  // around it right now.  FIXME! XXX GCW 2002-03-28

  model* m = const_cast<model*>(&mod);

  m->mean_cov_sub_model  .calculate_covariate_means(my_ped_data);
  m->var_cov_sub_model   .calculate_covariate_means(my_ped_data);
  m->susc_cov_sub_model  .calculate_covariate_means(my_ped_data);

  for(size_t trait_num = 0; trait_num < mean_cov_data.size(); ++trait_num)
  {
    double mean = mod.mean_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(mean_cov_data[trait_num], mean);
  }

  for(size_t trait_num = 0; trait_num < var_cov_data.size(); ++trait_num)
  {
    double mean = mod.var_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(var_cov_data[trait_num], mean);
  }

  for(size_t trait_num = 0; trait_num < susc_cov_data.size(); ++trait_num)
  {
    double mean = mod.susc_cov_sub_model.get_covariate_mean(trait_num);

    center_trait(susc_cov_data[trait_num], mean);
  }
}

//======================================================================
//                                                                     
//  calculate_transf_age_onset()
//                                                                     
//======================================================================

int
onset_member_calculator::calculate_transf_age_onset()
{
  // Copy over the ages

  transf_age_onsets = age_onsets;

  // Transform it.

  if(!mod.transf_sub_model.transform(transf_age_onsets))
  {
    return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     
//  calculate_transf_age_exam()
//                                                                     
//======================================================================

int
onset_member_calculator::calculate_transf_age_exam()
{
  // Copy over the ages

  transf_age_exams = age_exams;

  // Transform it.

  if(!mod.transf_sub_model.transform(transf_age_exams))
  {
    return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     
//  calculate_expected_means()
//                                                                     
//======================================================================

int
onset_member_calculator::calculate_expected_means()
{
  // Use temp vector to store the data

  vector<vector<double> >& mean_data = temp_calc_space;

  // Populate the vector with the non-polygenic values

  for(size_t i = 0; i < member_classes.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
      mean_data[genotype][i] = calc_exp_mean(i, genotype);
  }

  // Transform the means

  if(!mod.transf_sub_model.transform(mean_data[index_AA]))
  {
    return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  if(!mod.transf_sub_model.transform(mean_data[index_AB]))
  {
    return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  if(!mod.transf_sub_model.transform(mean_data[index_BB]))
  {
    return segreg_errors::MCC_FAILED_TRANSFORM;
  }

  // If the age is not polygenically affected, we copy the mean data into
  // the expected means and we're done

  if(mod.ons_sub_model.m_option() != onset_sub_model::m_A)
  {
    expected_means[0] = mean_data;

    return segreg_errors::EVAL_OK;
  }

  // Otherwise, we have to add our polygenic mean adjustment to it

  for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
  {
    expected_means[k] = mean_data;

    double poly_mean = mod.fpmm_sub_model.mean(k);

    for(size_t i = 0; i < affections.size(); ++i)
    {
      if(!is_member_valid(i)) continue;

      for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
        expected_means[k][genotype][i] += poly_mean;
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calc_exp_mean(...)                                                 =
//                                                                     =
//======================================================================
double
onset_member_calculator::calc_exp_mean(
        size_t abs_mem_ref,
        genotype_index genotype)
{
  double expected_mean = 0.0,
         tau           = 0.0;

  expected_mean = mod.mean_sub_model.parameter(genotype);

  for(size_t i = 0; i < mod.mean_cov_sub_model.covariates().size(); ++i)
  {
    expected_mean += mod.mean_cov_sub_model.covariates()[i].coefficient *
                     mean_cov_data[i][abs_mem_ref];

    if(mod.mean_cov_sub_model.covariates()[i].has_interaction)
    {
      tau = mod.mean_cov_sub_model.covariates()[i].i_taus[genotype];

      expected_mean += tau *
                       mod.mean_sub_model.parameter(genotype) *
                       mean_cov_data[i][abs_mem_ref];
    }
  }

  return expected_mean;
}

//======================================================================
//                                                                     
//  calculate_expected_suscs()
//                                                                     
//======================================================================

int
onset_member_calculator::calculate_expected_suscs()
{
  // Use temp vector to store the data

  vector<vector<double> >& susc_data = temp_calc_space;

  // Populate the vector with the non-polygenic values

  for(size_t i = 0; i < affections.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
      susc_data[genotype][i] = calc_exp_susc(i, genotype);
  }

  // If the susc is not polygenically affected, we copy the susc data into
  // the expected suscs and we're done

  if(mod.ons_sub_model.m_option() != onset_sub_model::m_S)
  {
    expected_susceptibilities[0] = susc_data;

    return segreg_errors::EVAL_OK;
  }

  // Otherwise, we have to add our polygenic component to it

  for(size_t k = 0; k < mod.fpmm_sub_model.max_pgt(); ++k)
  {
    expected_susceptibilities[k] = susc_data;

    double poly_susc = mod.fpmm_sub_model.mean(k);

    for(size_t i = 0; i < affections.size(); ++i)
    {
      if(!is_member_valid(i)) continue;

      for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
        expected_susceptibilities[k][genotype][i] += poly_susc;
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calc_exp_susc(...)                                                 =
//                                                                     =
//======================================================================
double
onset_member_calculator::calc_exp_susc(
        size_t abs_mem_ref,
        genotype_index genotype)
{
  double expected_susc = 0.0,
         tau           = 0.0;

  expected_susc = mod.susc_sub_model.parameter(genotype);

  for(size_t i = 0; i < mod.susc_cov_sub_model.covariates().size(); ++i)
  {
    expected_susc += mod.susc_cov_sub_model.covariates()[i].coefficient *
                     susc_cov_data[i][abs_mem_ref];

    if(mod.susc_cov_sub_model.covariates()[i].has_interaction)
    {
      tau = mod.susc_cov_sub_model.covariates()[i].i_taus[genotype];

      expected_susc += tau *
                       mod.susc_sub_model.parameter(genotype) *
                       susc_cov_data[i][abs_mem_ref];
    }
  }

  return expected_susc;
}

//======================================================================
//                                                                     
//  calculate_expected_alphas()
//                                                                     
//======================================================================

int
onset_member_calculator::calculate_expected_alphas()
{
  // Populate the expected alphas vector with the individual specific
  // values

  for(size_t i = 0; i < affections.size(); ++i)
  {
    if(!is_member_valid(i)) continue;

    for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
    {
      expected_alphas[genotype][i] = calc_exp_alpha(i, genotype);

      if(expected_alphas[genotype][i] < 0.0)
        return segreg_errors::MCC_VARIANCE_INVALID;
    }
  }

  return segreg_errors::EVAL_OK;
}

//======================================================================
//                                                                     =
//  calc_exp_alpha(...)                                                 =
//                                                                     =
//======================================================================
double
onset_member_calculator::calc_exp_alpha(
        size_t abs_mem_ref,
        genotype_index genotype)
{
  double expected_var = 0.0,
         tau           = 0.0;

  expected_var = mod.var_sub_model.parameter(genotype);

  for(size_t i = 0; i < mod.var_cov_sub_model.covariates().size(); ++i)
  {
    expected_var += mod.var_cov_sub_model.covariates()[i].coefficient *
                    var_cov_data[i][abs_mem_ref];

    if(mod.var_cov_sub_model.covariates()[i].has_interaction)
    {
      tau = mod.var_cov_sub_model.covariates()[i].i_taus[genotype];

      expected_var += tau *
                      mod.var_sub_model.parameter(genotype) *
                      var_cov_data[i][abs_mem_ref];
    }
  }

  double expected_alpha = M_PI / (sqrt(3.0) * sqrt(expected_var));

  return expected_alpha;
}

}}
