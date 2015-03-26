//============================================================================
// File:      transmission sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/1/01 - created.                               djb
//                                                          
// Notes:     transmission_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/transmission_sub_model.h"

namespace SAGE
{
namespace SEGREG
{

const std::string  TRANSMISSION_NAME     = "Transmission";
const double       TRANSM_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN();               // .5;                 
const double       TRANSM_HOMOG_GENERAL_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // .3;    
const double       TRANSM_HOMOG_MEND_AA_DEFAULT_VALUE = 1;
const double       TRANSM_HOMOG_MEND_AB_DEFAULT_VALUE = .5;
const double       TRANSM_HOMOG_MEND_BB_DEFAULT_VALUE = 0;
const double       TRANSM_AB_FREE_AA_DEFAULT_VALUE = 1;
const double       TRANSM_AB_FREE_BB_DEFAULT_VALUE = 0;
const bool         TRANSM_DEFAULT_FIXED = false;
const double       TRANSM_LB = 0;
const double       TRANSM_UB = 1;

//============================================================================
// IMPLEMENTATION:  transmission_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
// - Call the appropriate sub-model setting function for the 
//   specified sub-model option.  A return value of false means
//   that sub-model set failed.
//
bool
transmission_sub_model::set(sm_option opt, const model_input& tau_AA,
                            const model_input& tau_AB, const model_input& tau_BB, 
                            bool no_bounds, bool type_missing)
{
  if(isLinked())
  {
    return false;
  }

  // - Not specifying mean sub-model is a way of specifying commingling
  //   analysis, but certain constraints apply.  See also set functions
  //   in variance and genotype frequency sub-models.
  //
  if(type_missing)
  {
    my_errors << priority(error) << "transmission sub-block specified without "
              << "specifying a type (mean or suscept) sub-block.  Ignoring ..." << endl;
    return true;
  }
  
  // - Make sure there are no inputs for which fixed == true and
  //   value == QNAN.
  //
  // - Parser revised 8-28-01 at request of rce to dissallow situation
  //   where fixed == true and no value spec'ed for val except for
  //   a couple of special cases (not in this sub-model) so this
  //   call should always result in a return value of true.  -djb
  //
  if(! inputs_ok(tau_AA, tau_AB, tau_BB))
  {
    return false;
  }

  // Determine the bounds on transmission values.
  
  double  transm_lb = TRANSM_LB;
  double  transm_ub = TRANSM_UB;
  
  // - no bounds specification applies only to homog_general and general
  //   options.
  //
  if(opt == homog_general || opt == general)
  {
    if(no_bounds)
    {
      transm_lb = NEGATIVE_INF;
      transm_ub = POSITIVE_INF;
    }
  }
  else
  {
    if(no_bounds)
    {
      my_errors << priority(error) << "parameter, no_bounds, not relevant "
                << "for " << option_2_parameter(opt) << " option in transmission "
                << "sub-block.  Ignoring ..." << endl;
    }
  } 
  
  my_option = opt;

  bool return_value = false;
  switch(opt)
  {
    case no_trans:
      return_value = set_no_trans(tau_AA, tau_AB, tau_BB);
      break;
    case homog_no_trans:
      return_value = set_homog_no_trans(tau_AA, tau_AB, tau_BB);
      break;
    case homog_mendelian:
      return_value = set_homog_mendelian(tau_AA, tau_AB, tau_BB);
      break;
    case homog_general:
      return_value = set_homog_general(tau_AA, tau_AB, tau_BB, transm_lb, transm_ub);
      break;
    case general:
      return_value = set_general(tau_AA, tau_AB, tau_BB, transm_lb, transm_ub);
      break;
    case tau_ab_free:
      return_value = set_tau_ab_free(tau_AA, tau_AB, tau_BB, transm_lb, transm_ub);
      break;
    case mitochondrial:
      return_value = set_mitochondrial(tau_AA, tau_AB, tau_BB);
    default:
      ;       // Return value is false.  
  }
  
  return return_value;
}

bool
transmission_sub_model::set_no_trans
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB)    
{
  bool  return_value = false;

  // There are no parameters in the no_trans case.
  
  my_parameters.resize(0);

  // Do our setting based upon the frequency model
  
  assert(my_gf_ptr != 0);
  switch(my_gf_ptr->option())
  {
    case genotype_frequency_sub_model::hwe:

    case genotype_frequency_sub_model::NONE:
      return_value = set_no_trans_hwe(tau_AA, tau_AB, tau_BB);
      break;

    case genotype_frequency_sub_model::nhwe:
      return_value = set_no_trans_nhwe(tau_AA, tau_AB, tau_BB);
      break;

    default:
      assert(false);
  }
  
  return return_value;
}

// - If genotype frequency sub model option is hwe the transmission
//   sub model no_trans option is the same as homog_no_trans.
//
bool
transmission_sub_model::set_no_trans_hwe
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB)    
{
  // - SCR 190.
  //
  my_option = homog_no_trans;
  
  return set_homog_no_trans(tau_AA, tau_AB, tau_BB);
}

// - If genotype frequency sub model option is nhwe there are no
//   transmission probabilities (genotype probabilities conditional
//   on parent's genotypes same as population genotype probabilities).
//
bool
transmission_sub_model::set_no_trans_nhwe
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB)    
{
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);

  if(info)
  {
    my_errors << priority(error) << "Transmission probability specified in "
              << "transmission sub-block with no_trans option.  "
              << "Ignoring ..." << endl;
  }
  
  // Set all our taus to QNAN
  
  my_taus[index_AA] = QNAN;
  my_taus[index_AB] = QNAN;
  my_taus[index_BB] = QNAN;
  
  return true;
}


// - All taus equal frequency of allele A.
//
bool
transmission_sub_model::set_homog_no_trans
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB)
{
  // No parameters in the homog_no_trans case

  my_parameters.resize(0);

  // If there is any taus given, this indicates an error.
  
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);

  if(info)
  {
    my_errors << priority(error) << "Transmission probability specified in "
              << "transmission sub-block with no_trans or homog_no_trans option.  "
              << "Ignoring ..." << endl;
  }
  
  // Set all our taus = qA.
  
  copy_freq_a_into_taus();
  
  return true;
}

bool
transmission_sub_model::set_mitochondrial
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB)
{
  // No parameters in the mitochondrial case

  my_parameters.resize(0);

  // If there is any taus given, this indicates an error.
  
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);

  if(info)
  {
    my_errors << priority(error) << "Transmission probability specified in "
              << "transmission sub-block with mitochondrial option.  "
              << "Ignoring ..." << endl;
  }
  
  // Set all our taus = QNAN
  
  my_taus[index_AA] = QNAN;
  my_taus[index_AB] = QNAN;
  my_taus[index_BB] = QNAN;
  
  return true;
}


int
transmission_sub_model::finalizeConfiguration()
{
  if(my_option == no_trans       || 
     my_option == homog_no_trans ||
     my_option == mitochondrial)
  {
    my_parameters.resize(0);
  }
  else
  {
    finalize_parameters();
  }
  
  return 0;
}

void
transmission_sub_model::initialize_parameters
  (MAXFUN::Parameter::ParamTypeEnum AA_type,
   MAXFUN::Parameter::ParamTypeEnum AB_type,
   MAXFUN::Parameter::ParamTypeEnum BB_type,
   double                           lower_bound,
   double                           upper_bound
  )
{
  // Our model has three parameters, one for each of the taus.

  my_parameters.resize(3);

  // Define our initialization parameters based upon the 
  // arguments and the (already set) taus

  my_parameters[0] =
      MAXFUN::ParameterInput
        ( "TRANSMISSIONS", "transm prob_AA", AA_type,
          my_taus[index_AA], lower_bound, upper_bound);
  my_parameters[1] = 
      MAXFUN::ParameterInput
        ( "TRANSMISSIONS", "transm prob_AB", AB_type,
          my_taus[index_AB], lower_bound, upper_bound);
  my_parameters[2] =
      MAXFUN::ParameterInput
        ( "TRANSMISSIONS", "transm prob_BB", BB_type,
          my_taus[index_BB], lower_bound, upper_bound);
}

void
transmission_sub_model::finalize_parameters()
{
  my_parameters[index_AA].initial_estimate = my_taus[index_AA];
  my_parameters[index_AB].initial_estimate = my_taus[index_AB];
  my_parameters[index_BB].initial_estimate = my_taus[index_BB];
}

// - Taus AA, AB and BB are fixed at 1, .5 and 0 respectively.
//
bool
transmission_sub_model::set_homog_mendelian(const model_input& tau_AA, 
                                           const model_input& tau_AB, const model_input& tau_BB)
{
  // Specify our taus
  
  my_taus[index_AA] = TRANSM_HOMOG_MEND_AA_DEFAULT_VALUE;
  my_taus[index_AB] = TRANSM_HOMOG_MEND_AB_DEFAULT_VALUE;
  my_taus[index_BB] = TRANSM_HOMOG_MEND_BB_DEFAULT_VALUE;

  // Set up our initial my_parameters.

  initialize_parameters
      (MAXFUN::Parameter::FIXED,
       MAXFUN::Parameter::FIXED,
       MAXFUN::Parameter::FIXED,
       TRANSM_LB,
       TRANSM_UB);
  
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);
  
  if(info)
  {
    my_errors << priority(error) << "Transmission probability specified in "
              << "transmission sub-block with homog_mendelian option.  "
              << "Ignoring ..." << endl;
  }
  
  return true;
}


// - Tau AB is a function of taus AA and BB and qA.  User specified tau AB
//   is ignored.
//
bool
transmission_sub_model::set_homog_general
      (const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB, 
       double transm_lb, double transm_ub)
{ 
  // Specify our taus using default arguments.
  
  my_taus[index_AA] = TRANSM_HOMOG_GENERAL_DEFAULT_VALUE;
  my_taus[index_AB] = TRANSM_HOMOG_GENERAL_DEFAULT_VALUE;
  my_taus[index_BB] = TRANSM_HOMOG_GENERAL_DEFAULT_VALUE;

  // Initialize my_parameters with default arguments.
  
  initialize_parameters(
    MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
    MAXFUN::Parameter::DEPENDENT,
    MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
    transm_lb, transm_ub);
  
  // Determine our info and change any tau values that need changing due to it.

  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);

  // Initially assume our return is false.

  bool return_value = false;
  
  // Deal with the AB, which is ignored in homog_general
  
  if(info & AB)
  {
    my_errors << priority(error) << "Transmission probability AB specified in "
              << "transmission sub-block with homog_general option.  "
              << "Ignoring ..." << endl;
  }
  
  // Switch based upon just the AA and BB values
  
  switch(info & AA_BB)
  {
    case no_geno:
      // Note that this includes the just AB specified case.
      return_value = true;
      break;
      
    case AA:
      return_value = set_homog_general_one_value(tau_AA, true);
      break;
      
    case BB:
      return_value = set_homog_general_one_value(tau_BB, false);
      break;
      
      
    case AA_BB:
      return_value = set_homog_general_two_values(tau_AA, tau_BB);
      break;
      
    default:
      assert(false);
  }
  
  if(return_value)
  {
    // Since our values are good, we can calculate our tau_ab value now.
    
    calc_tau_ab();
    
    // Determine the state of tau_ab.  It is fixed if both tau_aa and tau_bb
    // are fixed, and dependent otherwise
    
    if(my_parameters[index_AA].initial_type == MAXFUN::Parameter::FIXED &&
       my_parameters[index_BB].initial_type == MAXFUN::Parameter::FIXED)
    {
      my_parameters[index_AB].initial_type = MAXFUN::Parameter::FIXED;
    }
    else
    {
      my_parameters[index_AB].initial_type = MAXFUN::Parameter::DEPENDENT;
    }

    // Our values are good, so update our maxfun initialization with our taus

    finalize_parameters();
  }
  
  return return_value;
}

bool
transmission_sub_model::set_general
      (const model_input& tau_AA, const model_input& tau_AB, 
       const model_input& tau_BB, double transm_lb, double transm_ub)
{ 
  // Start under the assumption that we're not good
  
  bool return_value = false;

  // Specify our taus using default arguments. 
  
  my_taus[index_AA] = TRANSM_DEFAULT_VALUE;
  my_taus[index_AB] = TRANSM_DEFAULT_VALUE;
  my_taus[index_BB] = TRANSM_DEFAULT_VALUE;

  // Initialize my_parameters with default arguments.
  
  initialize_parameters(
    MAXFUN::Parameter::INDEPENDENT,
    MAXFUN::Parameter::INDEPENDENT,
    MAXFUN::Parameter::INDEPENDENT,
    transm_lb, transm_ub);

  // Determine our information, and set our taus and their types
  
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);
  
  if(info == no_geno)
  {
    return_value = true;
  }
  else
  {
    if(info & AA)
    {
      return_value |= set_general_one_value(tau_AA, index_AA);
    }

    if(info & AB)
    {
      return_value |= set_general_one_value(tau_AB, index_AB);
    }

    if(info & BB)
    {
      return_value |= set_general_one_value(tau_BB, index_BB);
    }
  }
      
  if(return_value)
  {
    // Our values are good, so update our maxfun initialization with our taus

    finalize_parameters();
  }
  
  return return_value;
}

// - Taus AA and BB fixed at 1 and 0, respectively.  Tau AB
//   is settable.
//
bool
transmission_sub_model::set_tau_ab_free
      (const model_input& tau_AA, const model_input& tau_AB, 
       const model_input& tau_BB, double transm_lb, double transm_ub)
{ 
  // Specify our taus using default arguments. 
  
  my_taus[index_AA] = TRANSM_AB_FREE_AA_DEFAULT_VALUE;
  my_taus[index_AB] = TRANSM_DEFAULT_VALUE;
  my_taus[index_BB] = TRANSM_AB_FREE_BB_DEFAULT_VALUE;

  // Initialize my_parameters with default arguments.
  
  initialize_parameters(
    MAXFUN::Parameter::FIXED,
    MAXFUN::Parameter::INDEPENDENT,
    MAXFUN::Parameter::FIXED,
    transm_lb, transm_ub);

  // Determine our information, and set our taus and their types
  
  genotype_info  info = total_info(tau_AA, tau_AB, tau_BB);
  
  bool return_value = false;

  if(info == no_geno)
  {
    return_value = true;
  }
  else
  {
    // Ignore AA and BB values, if they're specified, but we still return true
    if(info & AA)
    {
      my_errors << priority(error) << "Transmission probability AA specified " 
                << "in transmission sub-block with tau_ab_free option.  Ignoring ..."
                << endl;

      return_value = true;
    }

    if(info & BB)
    {
      my_errors << priority(error) << "Transmission probability BB specified " 
                "in transmission sub-block with tau_ab_free option.  Ignoring ..." << endl;

      return_value = true;
    }

    // Set the AB value.  Note that it overrides any other return value.
    if(info & AB)
    {
      return_value = set_tau_ab_free_AB_value(tau_AB);
    }
  }
  
  if(return_value)
  {
    // Our values are good, so update our maxfun initialization with our taus

    finalize_parameters();
  }
  
  return return_value;
}

//
// ===============  Ancillary Functions
//

bool
transmission_sub_model::input_meets_constraints(const model_input& input) const
{
  return TRANSM_LB <= input.value && input.value <= TRANSM_UB;
}

// - Make sure there are no inputs for which fixed == true and 
//   value == QNAN.
//
bool
transmission_sub_model::inputs_ok(const model_input& tau_AA, const model_input& tau_AB, 
                                  const model_input& tau_BB)
{
  bool  return_value = false;
  if(! ( (tau_AA.fixed && SAGE::isnan(tau_AA.value)) ||
         (tau_AB.fixed && SAGE::isnan(tau_AB.value)) ||
         (tau_BB.fixed && SAGE::isnan(tau_BB.value))   ))
  {
    return_value = true;
  }
  else
  {
    my_errors << priority(critical) << "In transmission sub-block, a tau parameter specified "
              << "as 'fixed', but no 'val' specified.  Skipping analysis ..." << endl;
  }
  
  return return_value;
}

// - Calculate value of tau AB as a function of taus AA and BB and qA.
//
void  
transmission_sub_model::calc_tau_ab()
{
  assert(my_gf_ptr && my_parameters.size() == 3);

  // Get the frequency of allele A and it's complement.
  
  double  qA = my_gf_ptr->freq_A();
  double  compl_qA = 1.0 - qA;
  
  // - This equation should only be called for nwe option of genotype frequency sub-block,
  //   in which case, frequency of allele A is restricted so that assertions below are met.
  //   See user documentation transmission sub-block notes.
  //
  assert(qA != 0);
  assert(compl_qA != 0);

  double AA_adjust = qA * qA * my_taus[index_AA];
  double AB_adjust = compl_qA * compl_qA * my_taus[index_BB];

  double t_ab = (qA - AA_adjust) - AB_adjust;

  //lint -e{414} Cannot be 0.0
  t_ab /= (2.0 * qA * compl_qA);

  if(t_ab < 0.0) t_ab = 0.0;
  if(t_ab > 1.0) t_ab = 1.0;

  my_taus[index_AB] = t_ab;
}

// - Must only be called w. argument of tau AA or tau BB.
//
bool  
transmission_sub_model::set_homog_general_one_value(const model_input& mi, bool input_is_AA)
{
  bool  return_value = false;

  size_t  mi_index = input_is_AA ? 0 : 2;
  string  tau_name = input_is_AA ? "AA" : "BB";
  
  if(input_meets_constraints(mi))
  {
    my_default = false;

    // Set our tau and it's initial type.
    
    my_taus[mi_index] = mi.value;

    my_parameters[mi_index].initial_type =
        mi.fixed ? MAXFUN::Parameter::FIXED : 
                   MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;

    return_value = true;
  }
  else
  {
    if(mi.fixed)
    {
      my_errors << priority(critical) << "Transmission probability " << tau_name << " specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Skipping analysis ..." << endl;
    }
    else
    {
      my_errors << priority(error) << "Transmission probability " << tau_name << " specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Ignoring ..." << endl;    
      return_value = true;
    }
  }
  
  return return_value;
}


bool  
transmission_sub_model::set_homog_general_two_values
    (const model_input& mi_one, const model_input& mi_two)
{
  // - If either input does not meet constraints, situation is equivalent to a 
  //   one value set.
  //
  if(! input_meets_constraints(mi_one))
  {
    if(mi_one.fixed)
    {
      my_errors << priority(critical) << "Transmission probability AA specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Skipping analysis ..." << endl;
      return false;
    }
    else
    {
      my_errors << priority(error) << "Transmission probability AA specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Ignoring ..." << endl;
      return set_homog_general_one_value(mi_two, false);
    }
  }
  
  if(! input_meets_constraints(mi_two))
  {
    if(mi_two.fixed)
    {
      my_errors << priority(critical) << "Transmission probability BB specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Skipping analysis ..." << endl;
      return false;
    }
    else
    {
      my_errors << priority(error) << "Transmission probability BB specified in "
                << "transmission sub-block with homog_general option is not within "
                << "allowable range.  Ignoring ..." << endl;
      return set_homog_general_one_value(mi_one, true);
    }
  }
  
  // - The arguments are taus AA and BB and their values are w/i limits.
  //
  my_default = false;

  my_taus[index_AA] = mi_one.value;
  my_taus[index_BB] = mi_two.value;

  my_parameters[index_AA].initial_type = 
      mi_one.fixed ? MAXFUN::Parameter::FIXED :
                     MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
  my_parameters[index_BB].initial_type = 
      mi_two.fixed ? MAXFUN::Parameter::FIXED :
                     MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
  
  return true;
}
 
bool  
transmission_sub_model::set_general_one_value(const model_input& mi, genotype_index gindex)
{
  bool  return_value = false;

  if(input_meets_constraints(mi))
  {
    my_default = false;

    my_taus[gindex] = mi.value;

    my_parameters[gindex].initial_type =
        mi.fixed ? MAXFUN::Parameter::FIXED :
                   MAXFUN::Parameter::INDEPENDENT; 
    
    return_value = true;
  }
  else
  {
    if(mi.fixed == true)
    {
      my_errors << priority(critical) << "Transmission probability specified in "
                << "transmission sub-block with general option is not within "
                << "allowable range.  Skipping analysis ..." << endl;
    }
    else
    {
      my_errors << priority(error) << "Transmission probability specified in "
                << "transmission sub-block with general option is not within "
                << "allowable range.  Ignoring ..." << endl;
      return_value = true;
    }
  }
  
  return return_value;
}
 
bool  
transmission_sub_model::set_tau_ab_free_AB_value(const model_input& tau_AB)
{
  bool  return_value = false;

  if(input_meets_constraints(tau_AB))
  {
    my_default = false;

    my_taus[index_AB] = tau_AB.value;
    my_parameters[index_AB].initial_type =
        tau_AB.fixed ? MAXFUN::Parameter::FIXED :
                       MAXFUN::Parameter::INDEPENDENT; 

                       return_value = true;
  }
  else
  {
    if(tau_AB.fixed)
    {
      my_errors << priority(critical) << "Transmission probability AB specified in "
                << "transmission sub-block tau_ab_free option is not within "
                << "allowable range.  Skipping analysis ..." << endl;
      }
    else
    {
      my_errors << priority(error) << "Transmission probability AB specified in "
                << "transmission sub-block with tau_ab_free option is not within "
                << "allowable range.  Ignoring ..." << endl;
      return_value = true;
    }
  }
  
  return return_value;
}

//
// ==============  Synchronization w. Maxfun
//
int
transmission_sub_model::update()
{
  int  return_value = 1;

  switch(my_option)
  {
    case no_trans:
      return_value = synchronize_no_trans();
      break;
    case homog_no_trans:
      return_value = synchronize_homog_no_trans();
      break;
    case homog_mendelian:
      return_value = synchronize_homog_mendelian();
      break;
    case homog_general:
      return_value = synchronize_homog_general();
      break;
    case general:
      return_value = synchronize_general();
      break;
    case tau_ab_free:
      return_value = synchronize_tau_ab_free();
      break;
    case mitochondrial:
      return_value = synchronize_mitochondrial();
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_value;
}

// - Do nothing if genotype frequency sub-model option is nhwe; there
//   are no maxfun parameters and tau's were set to QNAN when sub-model
//   was set.
//
int
transmission_sub_model::synchronize_no_trans()
{
  assert(my_gf_ptr != 0);
  if(my_gf_ptr->option() == genotype_frequency_sub_model::hwe)
  {
    return synchronize_homog_no_trans();
  }
  else
  {
    return 0;
  }
}

// - No maxfun parameters of 'my own'.
//
int  
transmission_sub_model::synchronize_homog_no_trans()
{
  copy_freq_a_into_taus();
  return 0;
}

// - Synchronization not necessary since all parameters are fixed and internal
//   synchronization was done in set function.
//
int  
transmission_sub_model::synchronize_homog_mendelian()
{
  return 0;
}

// - Need to calculate tau AB.
//
int  
transmission_sub_model::synchronize_homog_general()
{
  // Get tau_AA and tau_BB from maxfun
  
  my_taus[index_AA] = getParam(index_AA);
  my_taus[index_BB] = getParam(index_BB);

  // Calculate tau_ab
  
  calc_tau_ab();
  
  // - Update maxfun w. newly calculated value of tau AB.
  //
  getParam(index_AB) = my_taus[index_AB];
  
  return 0;
}
 
int  
transmission_sub_model::synchronize_general()
{
  // Get taus from maxfun
  
  my_taus[index_AA] = getParam(index_AA);
  my_taus[index_AB] = getParam(index_AB);
  my_taus[index_BB] = getParam(index_BB);

  return 0;
}

int  
transmission_sub_model::synchronize_tau_ab_free()
{
  // Get tau_ab from maxfun
  
  my_taus[index_AB] = getParam(index_AB);

  return 0;
}

/// Synchronizes with maxfun in mitochondrial cases.
///
/// There are no maxfun parameters in this case.
int
transmission_sub_model::synchronize_mitochondrial()
{
    return 0;
}
 
}
}
