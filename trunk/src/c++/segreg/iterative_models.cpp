#include "segreg/iterative_models.h"
#include "segreg/freq_estimates.h"

// These defines are necessary to make life easier.  We undef them at the
// end of the file.

#define ONE          genotype_specific_sub_model::one
#define TWO          genotype_specific_sub_model::two
#define two_rec      genotype_specific_sub_model::two_rec
#define two_dom      genotype_specific_sub_model::two_dom
#define THREE        genotype_specific_sub_model::three
#define three_add    genotype_specific_sub_model::three_add
#define three_inc    genotype_specific_sub_model::three_inc
#define three_dec    genotype_specific_sub_model::three_dec

#define homog_no_trans  transmission_sub_model::homog_no_trans
#define homog_mendelian transmission_sub_model::homog_mendelian
#define no_trans        transmission_sub_model::no_trans
#define homog_general   transmission_sub_model::homog_general
#define tau_ab_free     transmission_sub_model::tau_ab_free
#define general         transmission_sub_model::general
#define mitochondrial   transmission_sub_model::mitochondrial

#define NONE            genotype_frequency_sub_model::NONE
#define hwe             genotype_frequency_sub_model::hwe
#define nhwe            genotype_frequency_sub_model::nhwe

namespace SAGE
{
namespace SEGREG
{

#include "im_helpers.cpp"

vector<primary_analysis_results> Iterative_Models::intermax; // due to JA for intermediate maximization

template<int M, transmission_sub_model::sm_option T> 
const primary_analysis_results&
Iterative_Models::get_model_results (mean_option m, bool quality)
{
  model_results& ms = my_iterative_models[model_classification(m,T)];

  // If we have already calculated this result, there's nothing to do.

  if(!ms.is_available())
  {
    // Create the model vector.

    model_vector v;

    get_starting_models<M,T>(v, m);

    // Indicate that we're starting the new model

    output_max_line(m, T);

    // If we have models, we can run them.  Otherwise, we have a problem.

    if(v.size())
    {
      // Transform the initial models to the current model.

      transform_starting_models<M,T>(v, m);

      run_models(m, T, ms, v, quality);
    }
    else
    {
      // We have no valid starting points.  We've failed to maximize here.

      produce_model_failure(m,T,ms);
    }

    // Indicate that we're finishing the model

    output_max_line(m, T, true);
  }

  return ms.result;
}

template<int M, transmission_sub_model::sm_option T>
void Iterative_Models::transform_starting_models
  (model_vector& v, mean_option m)
{ }

template<>
void Iterative_Models::get_starting_models<1,homog_no_trans>
  (model_vector& v, mean_option m)
{
  v.push_back(get_initial_model());
}

template<> 
void Iterative_Models::transform_starting_models<1,homog_no_trans>
    (model_vector& v, mean_option)
{
  if(use_target_mean_model(ONE))
  {
    get_type_model(v[0]) = get_type_model(my_target_model);
  }
}

///  Creates a set of two mean starting points from the initial and one mean
///  models.  This is functionality common to homog_no_trans and mendelian
///  models.
template <>
void Iterative_Models::get_starting_models<2,homog_no_trans> 
    (model_vector& v, mean_option m)
{
  // First, make sure the vector is empty

  v.clear();

  // We have two starting states, inital values and the one mean model. 
  // However, the one mean model might be faulty, so we only reserve enough for
  // two rather than resize.

  v.reserve(2);

  // The initial model is always assumed to be good.  This might not be the
  // case, with bad initial estimates, but there's no way to check that at
  // this point.

  v.push_back(get_initial_model());

  // The one mean model is also added, but only if it's valid

  const primary_analysis_results& one_mdl = get_model_results<1,homog_no_trans>(ONE, false);

  if(one_mdl.is_valid())
    v.push_back(one_mdl.get_final_model());
}

template <>
void Iterative_Models::get_starting_models<2,homog_mendelian> 
    (model_vector& v, mean_option m)
{
  get_starting_models<2,homog_no_trans>(v, m);
}

/// Creates the starting models for the 2 type, mitochondrial cases
///
/// This is identical to what's done for 2 types, homog_no_trans or mendelian
/// models
template <>
void Iterative_Models::get_starting_models<2,mitochondrial> 
    (model_vector& v, mean_option m)
{
  get_starting_models<2,homog_no_trans>(v, m);
}

///  Modifies a set of two mean, homog_no_trans starting states to valid
///  starting states.
///
///  From the starting models given, we have three variants for starting
///  means and genotype frequencies, as defined in the SEGREG formula
///  document, appendix B, section B
///
///  The starting conditions can vary from one starting state (user gave all
///  necessary information) to 6 starting states (no information).  This
///  varies based upon the status of the mean and frequency sub_models.
///
///  Number based on mean and frequency sub_models:
///
///  - if mean and freq model both complete:  Use those values as starting
///    points
///
///  - if mean complete, but freq is not: Use the three (or two) frequency
///    estimates available from the genotype frequency estimates as starting
///    points.
///
///  - if mean is incomplete, but freq is complete: Use the three estimates
///    for mean as given in Appendix B, section B.
///
///  - if neither is complete: Use the six (or 4) estimates as given in App. 
///    B, sec B. 

template <>
void Iterative_Models::transform_starting_models<2,homog_no_trans>
    (model_vector& v, mean_option m)
{

  // The initial models have to be updated to two models.  If the target's
  // initial estimates have been specified by the user, we use those.

  bool use_target_mean = use_target_mean_model(m);
  bool use_target_freq = use_target_freq_model(m);

  if(use_target_mean)
    update_means_to_target(v);

  if(use_target_freq)
    update_freqs_to_target(v);

  // If either of the mean or freq were not available (most of the time), we
  // must expand our options using create_two_mean_estimates

  if(!use_target_mean || !use_target_freq)
    create_two_mean_estimates(v, m);

  // Special Case:  When performing a three_dec mean model option, we must
  // reverse all the means for the two models.  This also requires
  // reversing the frequencies.  All other elements should not require
  // adjustment.

  if(get_type_model(my_target_model).option() == three_dec)
  {
    for(size_t i = 0; i < v.size(); ++i)
    {
      // Get the current estimates
      double bAA = get_type_model(v[i]).parameter(index_AA);
      double bBB = get_type_model(v[i]).parameter(index_BB);

      // Make certain they're ordered so bAA is larger than bBB

      if(bAA < bBB)
      {
        // Reverse all the estimates

        get_type_model(v[i]).set(m, bBB, QNAN, bAA);

        if(!use_target_freq)
        {
          double qA  = v[i].freq_sub_model.freq_A();

          v[i].freq_sub_model.set(hwe,
                (1.0 - qA), QNAN, QNAN, QNAN, model_input(0.0, true), false);
        }
      }
    }
  }
}

template <>
void Iterative_Models::transform_starting_models<2,homog_mendelian> 
    (model_vector& v, mean_option m)
{
  transform_starting_models<2,homog_no_trans>(v,m);

  // For all our models, change to mendelian inheritance
  for(size_t i = 0; i < v.size(); ++i)
  {
    v[i].transm_sub_model.set(homog_mendelian, QNAN, QNAN, QNAN, false, false);
  }
}

template <>
void Iterative_Models::transform_starting_models<2,mitochondrial> 
    (model_vector& v, mean_option m)
{
  transform_starting_models<2,homog_no_trans>(v,m);

  // For all our models, change to mendelian inheritance
  for(size_t i = 0; i < v.size(); ++i)
  {
    v[i].transm_sub_model.set(mitochondrial, QNAN, QNAN, QNAN, false, false);
  }
}

///  In the two mean, homog_general and tau_ab_free cases, we use the two,
///  homog_no_trans and homog_mendelian models as our starting points
///  (barring user defined start).  From these starting points, we create
///  several different models, which differ based on transmission
///  frequencies.  These models are transmission model specific, but the
///  rest of the function is general

template <>
void Iterative_Models::get_starting_models<2,homog_general> 
    (model_vector& v, mean_option m)
{
  v.clear();

  // Creating the model vector. There are two starting conditions.  The
  // first is the two, homog_no_trans, the second the two, homog_mendelian.

  // Adding the two, homog_no_trans:

  // If the homog_no_trans model is valid, we add it.  It is considered
  // valid if this is not a MLM model without residuals, *and* it
  // generates ok.

  if(my_target_model.get_model_class() != model_MLM    ||
     my_target_model.resid_sub_model.has_residuals()   )
  {
    const primary_analysis_results& two_ntr = two(m, homog_no_trans, false);

    if(two_ntr.is_valid()) v.push_back(two_ntr.get_final_model());
  }
     
  // If the model is MLM, with no residuals, the two, no_trans is not a
  // valid starting state.  In this case, we use all the starting states
  // we would have for the homog_no_trans.

  else
  {
    get_starting_models<2,homog_no_trans>(v,m);

    transform_starting_models<2,homog_no_trans>(v,m);
  }

  // Adding the two, homog_mendelian:

  const primary_analysis_results& two_mdl = two(m, homog_mendelian, false);

  if(two_mdl.is_valid()) v.push_back(two_mdl.get_final_model());

  // Update the mean and freq option, if the user has supplied them.

  if(use_target_mean_model(m))
    update_means_to_target(v);

  if(use_target_freq_model(m))
    update_freqs_to_target(v);
}

template <>
void Iterative_Models::get_starting_models<2,tau_ab_free> 
    (model_vector& v, mean_option m)
{
  get_starting_models<2,homog_general>(v,m);
}

template <>
void Iterative_Models::transform_starting_models<2,homog_general> 
    (model_vector& v, mean_option m)
{
  transform_models_to_homog_general(v);
}

template <>
void Iterative_Models::transform_starting_models<2,tau_ab_free> 
    (model_vector& v, mean_option m)
{
  transform_models_to_tau_ab_free(v);
}

///  In the two mean, general case, we use the two, homog_general and
///  tau_ab_free models as our starting points (barring user defined start).
///
template <>
void Iterative_Models::get_starting_models<2,general> 
    (model_vector& v, mean_option m)
{
  v.clear();

  if(use_target_mean_model(m) && use_target_freq_model(m))
  {
    v.resize(1, get_initial_model());
  }
  else
  {
    // Otherwise, we have to go to the two, homog_general and tau_ab_free.

    const primary_analysis_results& hgn = get_model_results<2,homog_general> (m, false);
    const primary_analysis_results& tab = get_model_results<2,tau_ab_free>   (m, false);

    if(hgn.is_valid()) v.push_back(hgn.get_final_model());
    if(tab.is_valid()) v.push_back(tab.get_final_model());
  }
}

template <>
void Iterative_Models::transform_starting_models<2,general> 
    (model_vector& v, mean_option m)
{
  if(use_target_mean_model(m))
    update_means_to_target(v);
  
  if(use_target_freq_model(m))
    update_freqs_to_target(v);

  // Now we check the transmission.  If the user gave it to us, we don't
  // need to construct them.

  transform_models_to_general(v);
}

template<>
void Iterative_Models::get_starting_models<3,homog_no_trans>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, homog_no_trans);
}

template<>
void Iterative_Models::get_starting_models<3,homog_mendelian>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, homog_mendelian);
}

template<>
void Iterative_Models::get_starting_models<3,mitochondrial>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, homog_mendelian);
}

template<>
void Iterative_Models::get_starting_models<3,homog_general>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, homog_general);

  model_vector vtmp;

  if(my_target_model.get_model_class() != model_MLM    ||
     my_target_model.resid_sub_model.has_residuals()   )
  {
    const primary_analysis_results& hnt = get_model_results<3,homog_no_trans>  (m, false);

    if(hnt.is_valid()) vtmp.push_back(hnt.get_final_model());
  }

  const primary_analysis_results& hmn = get_model_results<3,homog_mendelian> (m, false);

  if(hmn.is_valid()) vtmp.push_back(hmn.get_final_model());
  
  if(!vtmp.empty())
  {
    transform_models_to_homog_general(vtmp);

    v.reserve(v.size() + vtmp.size());

    v.insert(v.end(), vtmp.begin(), vtmp.end());
  }
}

template<>
void Iterative_Models::get_starting_models<3,tau_ab_free>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, tau_ab_free);

  const primary_analysis_results& hmn = get_model_results<3,homog_mendelian>  (m, false);

  if(hmn.is_valid())
  {
    model_vector vtmp(1, hmn.get_final_model());

    transform_models_to_tau_ab_free(vtmp);

    v.reserve(v.size() + vtmp.size());

    v.insert(v.end(), vtmp.begin(), vtmp.end());
  }
}

template<>
void Iterative_Models::get_starting_models<3,no_trans>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, no_trans);

  const primary_analysis_results& hnt = get_model_results<3,homog_no_trans>  (m, false);

  if(hnt.is_valid())
  {
    model_vector vtmp(1, hnt.get_final_model());

    vtmp[0].transm_sub_model.set(no_trans, QNAN, QNAN, QNAN, false, false);

    v.reserve(v.size() + vtmp.size());

    v.insert(v.end(), vtmp.begin(), vtmp.end());
  }
}

template<>
void Iterative_Models::get_starting_models<3,general>
  (model_vector& v, mean_option m)
{
  get_three_common_models(v, m, general);

  model_vector vtmp;

  if(my_target_model.get_model_class() != model_MLM    ||
     my_target_model.resid_sub_model.has_residuals()   )
  {
    const primary_analysis_results& ntr = get_model_results<3,no_trans>  (m, false);

    if(ntr.is_valid()) vtmp.push_back(ntr.get_final_model());
  }

  const primary_analysis_results& hgn = get_model_results<3,homog_general> (m, false);
  const primary_analysis_results& tab = get_model_results<3,tau_ab_free>   (m, false);

  if(hgn.is_valid()) vtmp.push_back(hgn.get_final_model());
  if(tab.is_valid()) vtmp.push_back(tab.get_final_model());
  
  if(!vtmp.empty())
  {
    transform_models_to_general(vtmp);

    v.reserve(v.size() + vtmp.size());

    v.insert(v.end(), vtmp.begin(), vtmp.end());
  }
}

template<>
void Iterative_Models::transform_starting_models<3,homog_no_trans>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,homog_mendelian>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,mitochondrial>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,homog_general>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,no_trans>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,tau_ab_free>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}

template<>
void Iterative_Models::transform_starting_models<3,general>
  (model_vector& v, mean_option m)
{
  transform_three_means(v,m);
}


void Iterative_Models::transform_three_means
  (model_vector& v, mean_option m) const
{
  if(use_target_mean_model(m))
    update_means_to_target(v);
  
  if(use_target_freq_model(m))
    update_freqs_to_target(v);
}

inline
void Iterative_Models::output_max_line
    (mean_option m, transm_option t, bool done) const
{
  if(m == ONE) t = homog_no_trans;

  out.messages() << "    Maximizing";

  size_t c = 0;

  if(output_mean)
  {
    string sing = singular(my_target_model);
    string plur = plural(my_target_model);

    string mean_string;

    switch(m)
    {
      case ONE       : mean_string = " one " + sing;              break;
      case TWO       : mean_string = " two " + plur;              break;
      case two_rec   : mean_string = " two " + plur + ", rec.";   break;
      case two_dom   : mean_string = " two " + plur + ", dom.";   break;
      case THREE     : mean_string = " three " + plur;            break;
      case three_add : mean_string = " three " + plur + ", add."; break;
      case three_inc : mean_string = " three " + plur + ", inc."; break;
      case three_dec : mean_string = " three " + plur + ", dec."; break;
    }

    out.messages() << mean_string;
    
    c += mean_string.size();
  }

  if(output_mean && output_transm)
  {
    out.messages() << ",";
    c+= 1;
  }
  
  if(output_transm)
  {
    string transm_string;

    switch(t)
    {
      case homog_no_trans  : transm_string = " homog. no transm.";      break;
      case homog_mendelian : transm_string = " mendelian transm.";      break;
      case no_trans        : transm_string = " no transm.";             break;
      case homog_general   : transm_string = " homog. general transm."; break;
      case tau_ab_free     : transm_string = " tau_ab free transm.";    break;
      case general         : transm_string = " general transm.";        break;
      case mitochondrial   : transm_string = " mitochondrial transm."; break;
    }

    out.messages() << transm_string;
    
    c += transm_string.size();
  }

  c = 48 - c;
  
  // There's an easier way to do this, but I can't remember *sigh*
  for(size_t i = 0; i < c; ++i)
    out.messages() << "."; 
  
  if(done) out.messages() << "...Done.";
  
  out.messages() << endl;
}

const primary_analysis_results& Iterative_Models::get_model
    (mean_option m, transm_option t, bool quality)
{
  switch(m)
  {
    case ONE       : return get_model_results<ONE,homog_no_trans>(ONE, quality);

    case TWO       :
    case two_dom   :
    case two_rec   : return two(m, t, quality);

    case THREE     :
    case three_add :
    case three_inc :
    case three_dec : return three(m, t, quality);

    default        : break;
  }      

  // SHOULD NEVER HAPPEN!!!
  out.errors() << priority(critical) << "Fatal Error in model "
               << "selection routines.  Please contact S.A.G.E. support." << endl;
  exit(1);

  // This statement is unreachable, and included only because the compiler
  // (KCC) gets annoyed when it isn't
  return get_model_results<ONE,homog_no_trans>(ONE, quality);
}

void Iterative_Models::clear()
{
  my_target_model  = model();
  my_initial_model = model();

  initial_model_built     = false;
  target_mean_complete    = false;
  target_transm_complete  = false;
  target_freq_complete    = false;

  my_iterative_models = std::map<model_classification, model_results>();
}

void Iterative_Models::create_initial_model(model& m)
{
  bool interrupt = false;

  if(output)
    out.messages() << "    Setting up initial model.................................." << flush;
  
  // Make the initial model into the target model

  my_initial_model = my_target_model;

  // The initial model is restricted to one type, homog_no_trans, one
  // variance, and no genotype_frequency. We now modify the model given to
  // fit these constraints.  We also fill in missing values that may occur.
  
  interrupt |= create_initial_type   (m);
  interrupt |= create_initial_var    (m);
  interrupt |= create_initial_transm (m);
  interrupt |= create_initial_freq   (m);
  interrupt |= create_initial_fpmm   (m);
  interrupt |= create_initial_onset  (m);

  if(interrupt)
    out.messages() << endl
                   << "    Setting up initial model.................................."
                   << flush;

  if(output)
    out.messages() << "...Done." << endl;
  
  initial_model_built = true;
}

bool Iterative_Models::create_initial_type(model& m)
{
  bool interrupt = false;

  mean_option opt = get_type_model(m).option();

  // If option is not one or there are missing means, we have to set initial
  // values

  if(opt != ONE || SAGE::isnan(get_type_model(m).parameter(index_AA)))
  {
    double yb = my_trait_sample.get_mean();
    
    // If a binary trait, must modify the mean appropriately
    if(is_model_binary(m))
    {
      if (yb < 1) {
                  yb = log(yb / (1.0 - yb));
                  }
         else {yb = 0.9;}
    }

    //lint -e(534) since we check later for failures
    get_type_model(m).set(ONE, yb,yb,yb);
    
    if(output && opt != ONE)
    {
      interrupt = true;

      out.messages() << endl << "      Setting to one mean"
                     << "........................................Done." << flush;
    }
  }

  return interrupt;
}

bool Iterative_Models::create_initial_var(model& m, double var)
{
  bool interrupt = false;

  // We only do the variance if we're not dealing with a discrete trait
  if(m.get_primary_trait_type() == pt_BINARY) return false;

  // We also delay evaluation of the variance for age of onset.  If the
  // variance has been given us (from create_initial_onset) we set it
  // anyway
  if(SAGE::isnan(var) && is_model_onset(m)) return false;

  // If the variance model is not 1 variance or the variances are
  // unavailable, we must calculate the variance

  mean_option opt = m.var_sub_model.option();

  if(opt != ONE || SAGE::isnan(m.var_sub_model.parameter(index_AA)))
  {
    // If we don't hav a variance, we have to calculate it.
    
    if(SAGE::isnan(var))
    {
      var = my_trait_sample.get_var();

      switch(my_initial_model.get_model_class())
      {
        case model_A       : var *= 3;  break;
        case model_D       : var *= 30; break;
        case model_FPMM    : var *= 20; break;  // Continuous, FPMM 

        case model_MLM     :
        case model_INVALID :
        default            : SAGE_internal_error(); break;  // Should never happen!
      }
    }

    //lint -e(534) since we check later for failures
    m.var_sub_model.set(ONE, var,var,var, false, false);
    
    if(opt != ONE && output)
    {
      interrupt = true;

      out.messages() << endl << "      Setting to one variance"
                     << "....................................Done." << flush;
    }
  }

  return interrupt;
}

bool Iterative_Models::create_initial_transm(model& m)
{
  bool interrupt = false;

  // Set transmission to homog_no_trans

  if(m.transm_sub_model.option() != homog_no_trans)
  {
    //lint -e(534) since we check later for failures
    m.transm_sub_model.set
        (homog_no_trans, QNAN, QNAN, QNAN, false, false);
    
    if(output)
    {
      interrupt = true;

      out.messages() << endl << "      Setting to homog_no_trans"
                     << "..................................Done." << flush;
    }
  }

  return interrupt;
}

bool Iterative_Models::create_initial_freq(model& m)
{
  // Set geno_freq to none

  if(m.freq_sub_model.option() != NONE)
  {
    // XXX Note that currently the assortative mating value alpha is fixed at 0.  This
    // is due to the equations not being quite right yet.

    //lint -e(534) since we check later for failures
    m.freq_sub_model.set
        (NONE, QNAN, QNAN, QNAN, QNAN, model_input(0., true), true);
  }

  return false;
}

bool Iterative_Models::create_initial_fpmm(model& m)
{
  // If we're not an FPMM model, there's nothing to do here.
  if(m.get_model_class() != model_FPMM) return false;

  // We also delay evaluation of the variance for age of onset
  if(is_model_onset(m)) return false;

  // We divide up the FPMM calculations into onset and non-onset varieties

  if(!m.fpmm_sub_model.is_complete())
  {
    double s = my_trait_sample.get_var();

    double factor = is_model_continuous(m) ? 10 : 2;

    m.fpmm_sub_model.set(factor * s, m.fpmm_sub_model.frequency(), m.fpmm_sub_model.loci());
  }

  return false;
}

bool Iterative_Models::create_initial_onset(model& m)
{
  bool interrupt = false;

  // If it's not an onset model, we can just leave.
  if(!is_model_onset(m)) return false;

  // We must calculate three things here: 
  //
  //   1. the variance of the age of onset (delayed from
  //      create_initial_var()),
  //
  //   2. the variance of the polygenic component (delayed from
  //      create_initial_fpmm()),
  //
  //   3. the mean/susc of the polygenic component
  //
  // To do this, we first calculate the mean and variance of the age of onset
  // and affection status.  Some of these values are precomputed in
  // my_trait_sample, depending on the type component.

  double aff_mean;
  double aff_var;
  double age_mean;
  double age_var;

  if(m.ons_sub_model.t_option() == onset_sub_model::t_S)
  {
    // The affection values come directly from my_trait_sample

    aff_mean = my_trait_sample.get_mean();
    aff_var  = my_trait_sample.get_var ();

    // age of onset values must be computed.  We use the age of onset unless
    // there are no values, in which case we use the age_at_exam.

    mean_split ms;
    
    ms(*my_ped_data.get_raw_data(), m.ons_sub_model.age_of_onset());

    if(!ms.get_n())
      ms(*my_ped_data.get_raw_data(), m.ons_sub_model.age_at_exam());

    age_mean = ms.get_mean();
    age_var  = ms.get_var ();
  }
  else // t_option == t_A 
  {
    // Affection values must be computed.

    mean_split ms;
    
    ms(*my_ped_data.get_raw_data(), m.get_primary_trait());

    aff_mean = ms.get_mean();
    aff_var  = ms.get_var ();

    // Age of onset values come directly from my_trait_sample
 
    age_mean = my_trait_sample.get_mean();
    age_var  = my_trait_sample.get_var ();
  }

  // Now that we have computed our means and variances, we can start setting
  // the parameters
  
  //   1. the variance of the age of onset

  bool mA = m.ons_sub_model.m_option() == onset_sub_model::m_A;
    
  double var = ((mA) ? 20 : 30) * age_var;
  
  interrupt |= create_initial_var(m, var);
  
  //   2. the variance of the polygenic component

  if(!m.fpmm_sub_model.is_complete())
  {
    bool mA = m.ons_sub_model.m_option() == onset_sub_model::m_A;
    double s;
    if (m.fpmm_sub_model.variance_fixed()) {
           s = m.fpmm_sub_model.variance();
     } 
     else 
     {
//       s = (mA)? 10*age_var : 2*aff_var; original code 
       s = (mA)? 2*age_var : 2*aff_var; // modified by JA
     }

     model_input minp(s,m.fpmm_sub_model.variance_fixed());
     m.fpmm_sub_model.set(minp, m.fpmm_sub_model.frequency(), m.fpmm_sub_model.loci());

// Original code here commented out by JA
//    bool mA = m.ons_sub_model.m_option() == onset_sub_model::m_A;
//
//    double s = (mA) ? 10 * age_var : 2 * aff_var;
//
//    m.fpmm_sub_model.set(s, m.fpmm_sub_model.frequency(), m.fpmm_sub_model.loci());

  }

  //   3. the mean/susc of the polygenic component

  if(m.ons_sub_model.t_option() == onset_sub_model::t_S)
  {
    // We must check the mean for validity.

    if(SAGE::isnan(m.mean_sub_model.parameter(index_AA)))
    {
      m.mean_sub_model.set(ONE, age_mean, age_mean, age_mean);
    }
  }
  else
  {
    // We must check the susceptibility for validity

    if(SAGE::isnan(m.susc_sub_model.parameter(index_AA)))
    {
      double susc = log(aff_mean / (1.0 - aff_mean));

      m.susc_sub_model.set(ONE, susc, susc, susc);
    }
  }

  return interrupt;
}

void Iterative_Models::transform_models_to_homog_general
    (model_vector& v) const
{
  
  // Updating the Transmission.  If the user has specified initial
  // values, we use those.  Otherwise, we must calculate some.

  if(use_target_transm_model(homog_general))
  {
    for(size_t i = 0; i < v.size(); ++i)
      v[i].transm_sub_model = my_target_model.transm_sub_model;
  }
  else
  {
    // We construct 4 different starting points: 
    //
    //     (0.9, 0.1), (0.9, 0.5),
    //     (0.5, 0.1), (0.5, 0.5)
    //
    // for each of the models we already have.  There are either 1 or 2
    // models in v at this point

    size_t initial_model_count = v.size();

    v.reserve(4 * initial_model_count);

    for(size_t i = 0; i < 3; ++i)
    {
      for(size_t j = 0; j < initial_model_count; ++j)
        v.push_back(v[j]);
    }
    
    for(size_t i = 0; i < initial_model_count; ++i)
    {
      v[i + 0 * initial_model_count].transm_sub_model.set(homog_general, 0.9, QNAN, 0.1, false, false);
      v[i + 1 * initial_model_count].transm_sub_model.set(homog_general, 0.9, QNAN, 0.5, false, false);
      v[i + 2 * initial_model_count].transm_sub_model.set(homog_general, 0.5, QNAN, 0.1, false, false);
      v[i + 3 * initial_model_count].transm_sub_model.set(homog_general, 0.5, QNAN, 0.5, false, false);
    }
  }
}

void Iterative_Models::transform_models_to_tau_ab_free
    (model_vector& v) const
{
  // Updating the Transmission.  If the user has specified initial
  // values, we use those.  Otherwise, we must calculate some.

  if(use_target_transm_model(tau_ab_free))
  {
    for(size_t i = 0; i < v.size(); ++i)
      v[i].transm_sub_model = my_target_model.transm_sub_model;
  }
  else
  {
    // We construct 2 different starting points: 
    //
    //     (qA), (0.5)
    //
    // for each of the models we already have.  There are either 1 or 2
    // models in v at this point

    size_t initial_model_count = v.size();

    v.reserve(2 * initial_model_count);

    for(size_t i = 0; i < initial_model_count; ++i)
      v.push_back(v[i]);
    
    for(size_t i = 0; i < initial_model_count; ++i)
    {
      double qA = v[i].freq_sub_model.freq_A();

      v[i + 0                  ].transm_sub_model.set(tau_ab_free, QNAN, qA,  QNAN, false, false);
      v[i + initial_model_count].transm_sub_model.set(tau_ab_free, QNAN, 0.5, QNAN, false, false);
    }
  }
}

void Iterative_Models::transform_models_to_general
    (model_vector& v) const
{
  // Now we check the transmission.  If the user gave it to us, we don't
  // need to construct them.

  if(use_target_transm_model(general))
  {
    for(size_t i = 0; i < v.size(); ++i)
      v[i].transm_sub_model = my_target_model.transm_sub_model;
  }
  else
  {
    // Set the transmission model to general, but don't change the taus.

    for(size_t i = 0; i < v.size(); ++i)
    {
      double tauAA = v[i].transm_sub_model.tau(index_AA);
      double tauAB = v[i].transm_sub_model.tau(index_AB);
      double tauBB = v[i].transm_sub_model.tau(index_BB);

      v[i].transm_sub_model.set(general, tauAA, tauAB, tauBB, false, false);
    }
  }
}

const primary_analysis_results& Iterative_Models::two
    (mean_option m, transm_option t, bool q)
{
  // We must choose which two model to build based on the transmission model

  switch(t)
  {
    case homog_no_trans   :
    case no_trans         : return get_model_results<2,homog_no_trans> (m, q);

    case homog_mendelian  : return get_model_results<2,homog_mendelian>(m, q);

    case homog_general    : return get_model_results<2,homog_general>  (m, q);

    case tau_ab_free      : return get_model_results<2,tau_ab_free>    (m, q);

    case general          : return get_model_results<2,general>        (m, q);
    
    case mitochondrial    : return get_model_results<2,mitochondrial>  (m, q);

    default               : // SHOULD NEVER HAPPEN
                            break;
  }
  SAGE_internal_error();

  // This statement is unreachable, and included only because the compiler
  // (KCC) gets annoyed when it isn't
  return get_model_results<2,homog_no_trans>(m, false);
}


const primary_analysis_results& Iterative_Models::three
    (mean_option m, transm_option t, bool q)
{
  // We must choose which two model to build based on the transmission model

  switch(t)
  {
    case homog_no_trans  : return get_model_results<3,homog_no_trans> (m,q);

    case homog_mendelian : return get_model_results<3,homog_mendelian>(m,q);

    case no_trans        : return get_model_results<3,no_trans>       (m,q);

    case homog_general   : return get_model_results<3,homog_general>  (m,q);

    case tau_ab_free     : return get_model_results<3,tau_ab_free>    (m,q);

    case general         : return get_model_results<3,general>        (m,q);

    case mitochondrial   : return get_model_results<3,mitochondrial>  (m,q);

    default              : break;
  }
  SAGE_internal_error();

  // This statement is unreachable, and included only because the compiler
  // (KCC) gets annoyed when it isn't
  return get_model_results<3,homog_no_trans>(m, false);
}

///  The three_basic option does three mean models for each of the many
///  different transmissions.  It is incomplete in that there may be
///  additional rules not explicitly specified in the three_basic that
///  are performed by the other three_* options.
///
///  The three_basic takes the appropriate two model and generates four new
///  starting states.  These starting states use any information specified by
///  the user (mean, transm, freq) or creates its own based upon the two
///  models.
void Iterative_Models::get_three_common_models 
  (model_vector& v, mean_option m, transm_option t)
{
  v.clear();

  const primary_analysis_results& two_results = two(TWO,t, false);

  if(two_results.is_valid())
  {
    v.push_back(two_results.get_final_model());

    // We must modify these here, since additional models that might be
    // added (in the general, homog_general, no_trans and tau_ab_cases) do
    // not need these modifications.

    if(use_target_transm_model(t))
    {
      v[0].transm_sub_model = my_target_model.transm_sub_model;
    }

    // If the mean is incomplete, we must expand it into four versions

    if(!use_target_mean_model(m))
    {
      create_three_mean_estimates(v, m);
    }

    // Adjust the variance & freqs for all models

    for(size_t i = 0; i < v.size(); ++i)
    {
      set_three_variance(v[i], !use_target_freq_model(m));

      // One special case:  For no_trans models, make sure the option is
      // no_trans, not homog_no_trans as in the two case.
      if(t == no_trans) 
        v[i].transm_sub_model.set(no_trans, QNAN, QNAN, QNAN, false, false);
    }
  }
}

/// Creates the initial estimates based upon the models given in v.

/// Based on whether the model is continuous or discrete, call the funtion
/// which creates initial estimates for mean, frequency, and variance for
/// two mean, homog_no_trans and homog_menelian models.

void Iterative_Models::create_two_mean_estimates
    (model_vector& v, mean_option two_opt)
{
// new adddition here to accomodate new starting values 
  if(is_model_continuous(v[0]))
  {
    genotype_frequency_initial_estimates::cont = true;
    create_continuous_two_mean_estimates(v, two_opt);
  }
  else // Discrete
  {
    genotype_frequency_initial_estimates::cont = false;
    create_discrete_two_mean_estimates(v, two_opt);
  }
}

/// Creates the initial estimates based upon the models given in v for continuous

/// Creates initial estimates for mean, frequency, and variance for
/// continous two mean, homog_no_trans and homog_menelian models.

void Iterative_Models::create_continuous_two_mean_estimates
    (model_vector& v, mean_option two_opt)
{
  
  if (Iterative_Models::intermax.size() != 0){
  one_mean_res = intermax[0].get_one_mean_results(); // due to JA for using
  }
 
  double beta_hat = one_mean_res;

  bool use_mean = !use_target_mean_model(two_opt);
  bool use_freq = !use_target_freq_model(two_opt);

  bool est1good = true;
  bool est2good = true;
  bool est3good = true;

  double bAA1,bAA2,bAA3,bBB1,bBB2,bBB3,q_A;

  //lint --e{534} <-- Suppress lint messages about ignored returns for this function

  // Create our frequency model.  
  
  genotype_frequency_initial_estimates
      freq_est(my_trait_sample, my_target_model.freq_sub_model);

// Needed to decide in advance which estimates to keep

  bAA1 = my_trait_sample.get_min();
  bAA2 = (my_trait_sample.get_sum() - my_trait_sample.get_max())/(my_trait_sample.get_n() - 1);
  bAA3 = my_trait_sample.get_y1bar();

  q_A = freq_est.get_model(0).freq_A();
  bBB1 = (beta_hat -q_A*q_A)/(1.0 -q_A*q_A);
  if (bBB1 < bAA1) est1good = false;

  q_A = freq_est.get_model(1).freq_A();
  bBB2 = (beta_hat -q_A*q_A)/(1.0 -q_A*q_A);
  if (bBB2 < bAA2) est2good = false;

  q_A = freq_est.get_model(2).freq_A();
  bBB3 = (beta_hat -q_A*q_A)/(1.0 -q_A*q_A);
  if (bBB3 < bAA3) est3good = false;

  size_t model_count = 0;

  if (est1good) model_count = model_count + 1;
  if (est2good) model_count = model_count + 1;
  if (est3good) model_count = model_count + 1;

  model_vector initial_models = v;

  // Clear v before creating our vector of models

  v.clear();

  // Insert our initial models model_count times

  size_t i;

  for(i = 0; i < model_count; ++i)
  {
    v.insert(v.begin(), initial_models.begin(), initial_models.end());
  }


  // We must now define the initial mean and freq estimates for these models, as

  // First set of estimates:

  size_t j =0;

  if (est1good) 
  {
    for(i = 0; i < initial_models.size(); ++i)
    {
      if(use_freq)
        v[i].freq_sub_model = freq_est.get_model(0);

      if(use_mean)
        get_type_model(v[i]).set(two_opt, bAA1, QNAN, bBB1);
    
      if(!use_mean && !use_freq)
        set_two_variance(v[i]);

   }
 }

 if (est1good) j = i;
 
  // Second set of estimates:

   if (est2good) {
    for(i = 0; i < initial_models.size(); ++i)
    {
      if(use_freq)
        v[i+j].freq_sub_model = freq_est.get_model(1);

      if(use_mean)
        get_type_model(v[j+i]).set(two_opt, bAA2, QNAN, bBB2);

      if(!use_mean && !use_freq)
        set_two_variance(v[j+i]);
    }
  }
   if (est2good) j = i+j;


  // Third set of estimates: 

  if(est3good) 
  {
    for(i = 0; i < initial_models.size(); ++i)
    {
      if(use_freq)
        v[j+i].freq_sub_model = freq_est.get_model(2);

      if(use_mean)
        get_type_model(v[j+i]).set(two_opt, bAA3, QNAN, bBB3);

      if(!use_mean && !use_freq)
        set_two_variance(v[j+i]);
    }
  }
  
}

void Iterative_Models::create_three_mean_estimates
    (model_vector& v, mean_option m) const
{
  // Resize our vector with the right number of initial estimates.  Generally,
  // this is 4, except with three_add models.  This is because the third and
  // fourth models use the AA and BB estmates from the two model directly,
  // so the AB is the same in both cases.
  v.resize(3 + (m != three_add), v[0]);
  
  // We must now set the initial estimates for our means in the 3 mean case. 
  // This is done according to the formulas in Appendix B.
  
  // Get the final estimates from the two value model.
  double bAA_hat = get_type_model(v[0]).parameter(index_AA);
  double bBB_hat = get_type_model(v[0]).parameter(index_BB);
  
  // For each option, we create an estimate for each of AA, AB, and BB
  // based upon the AA hat and BB hat above.
  //
  // NOTE:  the three_add model is a special case, as it does not need an
  //        AB to be specified.  Additionally, the third and fourth option
  //        which differ for other options, are identical under three_add,
  //        so it is removed.
  double AA_est, AB_est, BB_est;

  // First option
  // ------------
  //
  // This option takes makes the AA, AB and BB distributed such that the BB value
  // is the same as the two model, the AA and AB values averaged equals the AA 
  // from the two model, and all three values are spaced evenly.
  AA_est = (4.0 * bAA_hat - bBB_hat) / 3.0;
  AB_est = (m == three_add) ? QNAN : (2.0 * bAA_hat + bBB_hat) / 3.0;
  BB_est = bBB_hat;
  
  //lint -e(534) since we check later for failures
  get_type_model(v[0]).set(m, AA_est, AB_est, BB_est);
                           
  // Second option
  // ------------
  //
  // This option takes makes the AA, AB and BB distributed such that the AA value
  // is the same as the two model, the AB and BB values averaged equals the BB 
  // from the two model, and all three values are spaced evenly.
  AA_est = bAA_hat;
  AB_est = (m == three_add) ? QNAN : (2.0 * bBB_hat + bAA_hat) / 3.0; 
  BB_est = (4.0 * bBB_hat - bAA_hat) / 3.0;

  //lint -e(534) since we check later for failures
  get_type_model(v[1]).set(m, AA_est, AB_est, BB_est);

  // Third option
  // ------------
  //
  // This option lets AA and AB be as the AA from the two model, and BB be the same
  // as the BB from the two model
  AA_est = bAA_hat;
  AB_est = (m == three_add) ? QNAN : bAA_hat; 
  BB_est = bBB_hat;
  
  //lint -e(534) since we check later for failures
  get_type_model(v[2]).set(m, AA_est, AB_est, BB_est);

  // If we're doing a three_add model, we can end now
  if(m == three_add)
    return;

  // Fourth option
  // ------------
  //
  // This option lets AA and AB be as the AA from the two model, and BB be the same
  // as the BB from the two model
  AA_est = bAA_hat;
  AB_est = bBB_hat; 
  BB_est = bBB_hat;

  //lint -e(534) since we check later for failures
  get_type_model(v[3]).set(m, AA_est, AB_est, BB_est);
}



void Iterative_Models::run_models
    (mean_option m, transm_option t, model_results& ret, const model_vector& v, bool q)
{
  analyzer.set_quality(q);

  ret.result = analyzer.run_analysis(my_ped_data, v, false);

  if(ret.result.is_valid())
    ret.availability = model_results::good;
  else
    ret.availability = model_results::bad;
}

void Iterative_Models::produce_model_failure
    (mean_option m, transm_option t, model_results& mr)
{
  assert(!mr.result.is_valid());
  
  mr.availability = model_results::bad;

  out.errors() << priority(error) << "Unable to create valid initial "
                  "estimates for this model.  Check for outliers and/or reduce "
                  "the number of parameters to be estimated."
               << endl;
}


void Iterative_Models::set_two_variance(model& m)
{
  const genotype_frequency_sub_model& f = m.freq_sub_model;

  // Get the one mean variance

  double one_var = get_model_results<1,homog_no_trans>(ONE,false).get_final_model().
                     var_sub_model.parameter(index_AA);

  // Get everything else from the model in question

  double mean_AA = get_type_model(m).parameter(index_AA);
  double mean_AB = get_type_model(m).parameter(index_AB);
  double mean_BB = get_type_model(m).parameter(index_BB);

  double psi_AA = f.prob(index_AA);
  double psi_AB = f.prob(index_AB);
  double psi_BB = f.prob(index_BB);

  // Calculate the avg. mean (mu_bar)

  double avg_mean = mean_AA * psi_AA +
                    mean_AB * psi_AB +
                    mean_BB * psi_BB;

  // Calculate the new variance.

  double var = one_var - psi_AA * (mean_AA - avg_mean) * (mean_AA - avg_mean)
                       - psi_AB * (mean_AB - avg_mean) * (mean_AB - avg_mean)
                       - psi_BB * (mean_BB - avg_mean) * (mean_BB - avg_mean);

  // This is unlikely to ever happen, but just in case, let's
  // return something that's not at the bound.
  if(var < 1.e-6) var = 1.e-5;

  m.var_sub_model.set(ONE, var, var, var, false, false);
}

///  Sets the three mean frequencies and variances based on the two mean
///  equivalents such that:
///
///  -#  The mean of the two and three mean cases are the same (by adjusting
///      frequencies.
///  -#  The variance of the three mean cases is lower than the variance of
///      the two mean case.

void Iterative_Models::set_three_variance(model& m, bool use_freq)
{
  // Get the two mean equivalent model

  const model& tm = two(TWO, m.transm_sub_model.option(), false).get_final_model();

  double mean_A = get_type_model(tm).parameter(index_AA);
  double mean_B = get_type_model(tm).parameter(index_BB);

  double psi_A  = tm.freq_sub_model.prob(index_AA) + tm.freq_sub_model.prob(index_AB);
  double psi_B  = tm.freq_sub_model.prob(index_BB);

  double var2   = tm.var_sub_model.parameter(index_AA);

  double qA     = tm.freq_sub_model.freq_A();

  double mu_hat = mean_A * psi_A + mean_B * psi_B; 

  // Get means from the model in question

  double mean_AA = get_type_model(m).parameter(index_AA);
  double mean_AB = get_type_model(m).parameter(index_AB);
  double mean_BB = get_type_model(m).parameter(index_BB);

  double qA1;

  // Calculate qA1.  To do this, we calculate two values by quadratic
  // formula.  Often, however, the A is 0, so it's a simple equation to
  // solve

  double qA2;

  double A = mean_AA - 2.0 * mean_AB + mean_BB;

  // If A is close to 0, we can calculate qA1 directly
  if(abs(A) < 1.e-10)
  {
    qA1 = (mu_hat - mean_BB) / 2.0 / (mean_AB - mean_BB);
    qA2 = -1;
  }
  else
  {
    double B = 2.0 * (mean_AB - mean_BB);
    double C = mean_BB - mu_hat;

    double sqrt_term = B * B - 4 * A * C;

    sqrt_term = (sqrt_term < 0) ? 0 : sqrt(sqrt_term);

    qA1 = (-1.0 * B + sqrt_term) / 2.0 / A;
    qA2 = (-1.0 * B - sqrt_term) / 2.0 / A;
  }

  // Make qA1 the better option.  We'll ignore qA2 from here.  Better is
  // defined as the one closer to the two mean qA

  bool qA1_valid = 0.0 <= qA1  && qA1 <= 1.0;
  bool qA2_valid = 0.0 <= qA2  && qA2 <= 1.0;

  if(qA2_valid)
  {
    if(qA1_valid)
    {
      double dist1 = abs(qA1 - qA);
      double dist2 = abs(qA2 - qA);

      qA1 = (dist1 < dist2) ? qA1 : qA2;
    }
    else
      qA1 = qA2;
  }

  if(qA1   < 0.001) qA1 = 0.001;
  if(0.999 < qA1)   qA1 = 0.999;

  // Set the freq model.  Note that we determine here if hwe is used for
  // the three means or not.

  freq_option fo = my_target_model.freq_sub_model.option();

  m.freq_sub_model.set(fo, qA1, QNAN, QNAN, QNAN, model_input(0.0, true), false);

  // Calculate our three mean psi's

  double psi_AA = qA1 * qA1;
  double psi_AB = 2.0 * qA1 * (1.0 - qA1);
  double psi_BB = (1.0 - qA1) * (1.0 - qA1);

  // Calculate the variance adjustment, based on the two and three mean
  // frequencies and means.

  double v_adj = abs(
                     psi_A  * (mean_A  - mu_hat) * (mean_A  - mu_hat)
                   + psi_B  * (mean_B  - mu_hat) * (mean_B  - mu_hat)
                   - psi_AA * (mean_AA - mu_hat) * (mean_AA - mu_hat)
                   - psi_AB * (mean_AB - mu_hat) * (mean_AB - mu_hat)
                   - psi_BB * (mean_BB - mu_hat) * (mean_BB - mu_hat) 
                 );

  double var3 = var2 - v_adj;

  if(var3 < 1.e-6) var3 = 1.e-5;
                     
  m.var_sub_model.set(ONE, var3, var3, var3, false, false);

}


void Iterative_Models::create_discrete_two_mean_estimates
    (model_vector& v, mean_option two_opt)
{

   
// use_mean and use_freq are false if all genotype frequencies are known

  one_susc_res = intermax[0].get_one_susc_results(); // due to JA for using

  double beta_hat = one_susc_res; // Appendix B Eq.2 A.

  bool use_mean = !use_target_mean_model(two_opt);
  bool use_freq = !use_target_freq_model(two_opt);

  model_vector initial_models = v;

  genotype_frequency_initial_estimates
      freq_est(my_trait_sample, my_target_model.freq_sub_model);

  size_t model_count = 18;

  v.clear();

  // Insert our initial models model_count times

  for(size_t i = 0; i < model_count; ++i)
  {
    v.insert(v.begin(), initial_models.begin(), initial_models.end());
  }

   
  double f = std::numeric_limits<double>::quiet_NaN();
  double n = my_trait_sample.get_n();

  for (size_t k = 0; k < 5; ++k){
  for (size_t j = 0; j < 3; ++j){
  for (size_t i = 0; i < initial_models.size(); ++i){
    size_t index = (6*k + 2*j + i);
    if (j == 0) f = 1.0/n;
    if (j == 1) f = 0.2;
    if (j == 2) f = 0.4;
    if(use_freq) v[index].freq_sub_model = freq_est.get_model(k);
    double q_A = v[index].freq_sub_model.freq_A();
    if(use_mean)
    {
      double bAA = log(f/(1.0 -f));
      double bBB = (beta_hat -(q_A*q_A)*bAA)/(1.0 -q_A*q_A);
      get_type_model(v[index]).set(two_opt, bAA, QNAN, bBB);
    }
   }
  }
 }

  return;

}
}
}

#undef ONE
#undef two          
#undef two_rec      
#undef two_dom      
#undef three        
#undef three_add    
#undef three_inc    
#undef three_dec    

#undef homog_no_trans
#undef homog_mendelian
#undef no_trans
#undef homog_general
#undef tau_ab_free
#undef general

#undef NONE
#undef hwe
#undef nhwe
