//============================================================================
// File:        segreg_calculator.cpp
//
//
// Author:      Kai He
//
// History:     Aug 6, 2001
//
// Copyright (c) 2001 R. C. Elston
// All rights reserved
//============================================================================

/** @file
 *  This file has calculation methods for pedigree prob or likelihood and 
 *  individual prob or likelihood for different models such as calculate_A, 
 *  calculate_D3, and calculate_FMPP, and ascertainment calculation.
 *  The functions instantiates class model calculation objects and cooresponding
 *  function got called.
 *
 *  It creates pedigree section based on valid pedigree data, set start point,
 *  add each member into pedigree section and build it, then create regressive
 *  peeler or FPMM peeler based on model_class, set SL_calculator object for
 *  the peeler, and then call peeler functions: anterior and posterior, finally
 *  calculate prob.
 */

#include <functional>
#include <boost/mem_fn.hpp>
#include "segreg/segreg_calculator.h"
#include "segreg/segreg_utilities.h"
 
using namespace std;
namespace SAGE {
namespace SEGREG {

/// The PEN_CUTOFF is the lower bound on reported penetrance values. 
/// Penetrance values that are less than this cutoff are not stored.

static const double PEN_CUTOFF=1e-5;

vector<double> segreg_calculator::reg_ant_vec; // due to JA for likelihood rescaling
vector<double> segreg_calculator::reg_pos_vec; // see above
vector<double> segreg_calculator::mlm_ant_vec; // due to JA for likelihood rescaling
vector<double> segreg_calculator::mlm_pos_vec; // see above
vector<double> segreg_calculator::fpmm_ant_vec; // due to JA for likelihood rescaling
vector<double> segreg_calculator::fpmm_pos_vec; // see above

/// reduce_normal_set takes three double values, (AA, AB and BB), which are
/// assumed to sum to 1.0, removes all values less than the cutoff, then
/// re-normalizes, and returns the values.

inline void
reduce_normal_set(double& AA, double& AB, double& BB, double cutoff)
{
  if(AA < cutoff) AA = 0.0;
  if(AB < cutoff) AB = 0.0;
  if(BB < cutoff) BB = 0.0;

  double sum = AA + AB + BB;

  AA /= sum;
  AB /= sum;
  BB /= sum;
}

/// renormalize takes three double values, (AA, AB and BB), normalizes to
/// 1.0, removes all values less than the cutoff, then re-normalizes, and
/// returns the values.

inline void
renormalize(double& AA, double& AB, double& BB, double cutoff = PEN_CUTOFF)
{
  double sum = AA + AB + BB;

  AA /= sum;
  AB /= sum;
  BB /= sum;

  reduce_normal_set(AA, AB, BB, cutoff);

}

/// renormalize takes three log_double values, (AA, AB and BB), normalizes
/// to 1.0, removes all values less than the cutoff, then re-normalizes, and
/// returns the double equivalents

void renormalize(double&     AA, double&     AB, double&     BB,
                 log_double lAA, log_double lAB, log_double lBB,
                 double cutoff = PEN_CUTOFF)
{
  log_double sum(lAA + lAB + lBB);

  AA = (lAA / sum).get_double();
  AB = (lAB / sum).get_double();
  BB = (lBB / sum).get_double();

  reduce_normal_set(AA, AB, BB, cutoff);

}

static const long trunc_val = 1L<<31;

//----------------------------------------------------------------------------
//
//	evaluate(...)
//
//----------------------------------------------------------------------------

double segreg_calculator::evaluate(parameter_vector& v)
{
  // If we've already evaluated this, we just clear and return the
  // last_likelihood

  if(finite(last_likelihood))
  {
    double d = last_likelihood;

    last_likelihood = QNAN;

    return d;
  }

  // This should never happen!!

  cout << "ERROR: Possible bug in SEGREG.  Please Contact SAGE Immediately." << endl;
  exit(1);

  return internal_evaluate(v);
}
//----------------------------------------------------------------------------
//
//      internal_evaluate(...)
//
//---------------------------------------------------------------------------- 

#if WINDOWS || SUN || __INTEL
#define trunc(x) (x)
#endif

double segreg_calculator::internal_evaluate(parameter_vector& v)
{

#if 0
  static double lastvals[3];

//  if(nfe % 5 == 0)
    for(size_t i = 0; i < 3; ++i)
    {
      cout << i <<"\t" << doub2str(v[i],35) << ' '
           << doub2str(v[i]-lastvals[i],35) << endl;

      lastvals[i] = v[i];
    }

#endif


  double result = calculate().get_log();

  // We round off this number to 31 binary digits (approximately 10 digits)
  // past the decimal.  This is to make SEGREG somewhat more deterministic
  result = trunc(result * trunc_val) / trunc_val;

#if 0
  static double bestvalue = NEGATIVE_INF;

  if(nfe % 1 == 0)
    cout << "Result(" << nfe << ")\t=   " << doub2str(result, 35) 
         << ' ' << doub2str(bestvalue - result,35) << endl;

  if(result > bestvalue) bestvalue = result;

#endif

  ++nfe;

  return result;
}

//----------------------------------------------------------------------------
//
//      update_bounds(...)
//
//---------------------------------------------------------------------------- 

int segreg_calculator::update_bounds(parameter_vector& v)
{
   int err_code = 0;
   int asc_err_code = 0;

   model_class  mc = md.get_model_class();

   switch(mc)
   {
     case model_A    :
     case model_D    : err_code       = my_like_elts->update();
                       if(my_asc_like_elts)
                         asc_err_code = my_asc_like_elts->update();
                       break;
     case model_FPMM : err_code = fpmmsl->update();
                       if(asc_fpmmsl)
                         asc_err_code = asc_fpmmsl->update();

                       break;
     case model_MLM  : err_code = bmc->update();

                       if(!err_code && 
                          !my_mlm_corr_verifier->current_correlation_estimates_are_valid())
                         err_code = segreg_errors::MLM_CORR_BOUNDARY;

                       if(abmc)
                         asc_err_code = abmc->update();
                       break;
     default         : break;
   }

   // If there's an error, return it

   if(err_code)     return err_code;
   if(asc_err_code) return asc_err_code;

   // Ok, no errors so far.  Let's pull the wool over Maxfun's eyes

   last_likelihood = internal_evaluate(v);

   if(!finite(last_likelihood))
     return segreg_errors::BAD_LIKELIHOOD;

   return segreg_errors::EVAL_OK;
}

//----------------------------------------------------------------------------
//      
//      calculator(...)
//      
//      Return the sum of the subpedigree likelihoods for regressive model
//----------------------------------------------------------------------------

log_double segreg_calculator::calculate()
{  

  log_double prob(1.0);

  log_double continuous_penalty(exp(1.0));

  // Turned off until fix for two -> three mean problem resolved 
  //continuous_penalty.pow(calculate_continuous_penalty());

  //if(!finite(continuous_penalty.get_double()))
  
  continuous_penalty = 1.0;

  log_double prevalence_penalty(exp(1.0));

  // calculate prevalence penalty only for traits that aren't continuous
  // (binary and onset traits)

  primary_type pt = md.get_primary_trait_type();

  if(pt != pt_CONTINUOUS)
    prevalence_penalty = calculate_prevalence_penalty();

//  IF the prevalence is causing it to be infinities or QNaN, that should be
//  a bad likelihood.
//  if(!finite(prevalence_penalty.get_double()))
//    prevalence_penalty = 1.0;

  model_class  mc = md.get_model_class();

  // Calculate the unpenalized likelihood

  prob = calc_likelihood();

  // Add the penalizations

  prob *= continuous_penalty;

  if(mc == model_MLM || (mc == model_FPMM && pt != pt_CONTINUOUS))
    prob *= prevalence_penalty;

  return prob;
}

void segreg_calculator::calc_connected_regressive()
{
  for(PedigreeDataSet::SubpedigreeIterator
          subped = my_ped_data.get_subpedigree_begin();
          subped != my_ped_data.get_subpedigree_end(); ++subped)
  {
    // Get the general pedigree likelihood
    log_double like = calc_connected_regressive(*subped, false);

    // If ascertainment is turned on, we must adjust the pedigree likelihood
    if(using_ascertainment())
      like /= calc_connected_regressive(*subped, true);

    if(!accumulate_likelihood(like)) return;
  }
}

void segreg_calculator::calc_connected_FPMM()
{
  for(PedigreeDataSet::SubpedigreeIterator
          subped = my_ped_data.get_subpedigree_begin();
          subped != my_ped_data.get_subpedigree_end(); ++subped)
  {
    // Get the general pedigree likelihood

    log_double like = calc_connected_FPMM(*subped, false);

    // If ascertainment is turned on, we must adjust the pedigree likelihood
    if(using_ascertainment())
      like /= calc_connected_FPMM(*subped, true);

    if(!accumulate_likelihood(like)) return;
  }
}

void segreg_calculator::calc_connected_MLM()
{
  for(PedigreeDataSet::SubpedigreeIterator
          subped = my_ped_data.get_subpedigree_begin();
          subped != my_ped_data.get_subpedigree_end(); ++subped)
  {
    // Get the general pedigree likelihood
    log_double like = calc_connected_MLM(*subped, false);

    // If ascertainment is turned on, we must adjust the pedigree likelihood
    if(using_ascertainment())
      like /= calc_connected_MLM(*subped, true);

    if(!accumulate_likelihood(like)) return;
  }
}

log_double segreg_calculator::calc_connected_regressive
    (const FPED::Subpedigree& ps, bool asc)
{
  // Create and set up the regressive_peeler

  regressive_peeler plr(ps,  *((asc) ? my_asc_like_elts : my_like_elts));

  // Initially, the likelihood is 0

  log_double like(0.0);

  double ant_norm = 0; // due to JA for likelihood rescaling
  double pos_norm = 0; // due to JA for likelihood rescaling
  double rescale_like = 0;// see above

  segreg_calculator::reg_ant_vec.clear(); // due to JA for likelihod rescaling
  segreg_calculator::reg_pos_vec.clear(); // due to JA for likelihood rescaling

  // Doesn't matter which individual, so just use the first one


  const member_type& indi = ps.member_index(0);

  const TypeDescription& tdesc = my_like_elts->get_type_description();

  for(TypeDescription::StateIterator state = tdesc.begin();
      state != tdesc.end(); ++state)
  {
    like += calc_subped_given_member_regressive(&plr, indi, *state);
  }
  
 for (unsigned i = 0; i != segreg_calculator::reg_ant_vec.size(); i++) { // defining the rescaling constants
    ant_norm += segreg_calculator::reg_ant_vec[i];
     pos_norm += segreg_calculator::reg_pos_vec[i];
 }

 for (unsigned i = 0; i != segreg_calculator::reg_ant_vec.size(); i++) { // now rescale
    segreg_calculator::reg_ant_vec[i] = segreg_calculator::reg_ant_vec[i]/ant_norm;
     segreg_calculator::reg_pos_vec[i] = segreg_calculator::reg_pos_vec[i]/pos_norm;
 }

 for (unsigned i = 0; i != segreg_calculator::reg_ant_vec.size(); i++) { // the likelihood after rescaling 
    rescale_like += segreg_calculator::reg_ant_vec[i]*segreg_calculator::reg_pos_vec[i]; 
 }

 log_double new_log_like(0.0);
 new_log_like = ant_norm*pos_norm*rescale_like; // due to JA
  return new_log_like;

//  return like;
}

log_double segreg_calculator::calc_connected_FPMM
    (const FPED::Subpedigree& ps, bool asc)
{
  // Create and set up the FPMM_peeler

  FPMM_peeler* plr = new FPMM_peeler(ps, md);

  if(!asc)
  {
    plr->set_SL(fpmmsl);
    fpmmsl->set_peeler(plr);
  }
  else
  {
    plr->set_SL(asc_fpmmsl);
    asc_fpmmsl->set_peeler(plr);
  }

  // Initially, the likelihood is 0

  log_double like(0.0);

  double ant_norm = 0; // due to JA for likelihood rescaling
  double pos_norm = 0; // due to JA for likelihood rescaling
  double rescale_like = 0;// see above

  segreg_calculator::fpmm_ant_vec.clear(); // due to JA for likelihod rescaling
  segreg_calculator::fpmm_pos_vec.clear(); // due to JA for likelihood rescaling

  // Doesn't matter which individual, so just use the first one
  const member_type& indi = ps.member_index(0);

  genetic_info g;

  for(g.genotype = index_AA; g.genotype != index_INVALID; ++g.genotype)
  {
    for(g.polygenotype = 0; g.polygenotype < md.fpmm_sub_model.max_pgt(); ++g.polygenotype)
    {
      like += calc_subped_given_member_FPMM(plr, indi, g);
    }
  }

 for (unsigned i = 0; i != segreg_calculator::fpmm_ant_vec.size(); i++) { // defining the rescaling constants
    ant_norm += segreg_calculator::fpmm_ant_vec[i];
     pos_norm += segreg_calculator::fpmm_pos_vec[i];
 }

 for (unsigned i = 0; i != segreg_calculator::fpmm_ant_vec.size(); i++) { // now rescale
    segreg_calculator::fpmm_ant_vec[i] = segreg_calculator::fpmm_ant_vec[i]/ant_norm;
     segreg_calculator::fpmm_pos_vec[i] = segreg_calculator::fpmm_pos_vec[i]/pos_norm;
 }

 for (unsigned i = 0; i != fpmm_ant_vec.size(); i++) { // the likelihood after rescaling
    rescale_like += segreg_calculator::fpmm_ant_vec[i]*segreg_calculator::segreg_calculator::fpmm_pos_vec[i]; 
 }

  delete plr;
  

 log_double new_log_like(0.0);
 new_log_like = ant_norm*pos_norm*rescale_like; // due to JA
 
  return new_log_like;

//  return like;
}

log_double segreg_calculator::calc_connected_MLM
    (const FPED::Subpedigree& ps, bool asc)
{
  // Create and set up the mlm_peeler
  
  MlmLikelihoodElements lelt(*my_ped_data.get_raw_data(),md,asc);

  if(!asc)
    mpl = new mlm_peeler(ps,bmc,*my_ped_data.get_raw_data(),md,lelt, asc);
  else
    mpl = new mlm_peeler(ps,abmc,*my_ped_data.get_raw_data(),md,lelt,asc);

  // Initially, the likelihood is 0

  log_double like(0.0);

  // Initially, the likelihood is 0

  like = 0.0;

  double ant_norm = 0; // due to JA for likelihood rescaling
  double pos_norm = 0; // due to JA for likelihood rescaling
  double rescale_like = 0;// see above

  segreg_calculator::mlm_ant_vec.clear(); // due to JA for likelihod rescaling
  segreg_calculator::mlm_pos_vec.clear(); // due to JA for likelihood rescaling

  // Doesn't matter which individual, so just use the first one

  const member_type& indi = ps.member_index(0);

  for(TypeDescription::StateIterator state = lelt.get_type_description().begin();
      state != lelt.get_type_description().end(); ++state)
  {
    like += calc_subped_given_member_MLM(mpl, indi, *state);
  }

 for (unsigned i = 0; i != segreg_calculator::mlm_ant_vec.size(); i++) { // defining the rescaling constants
    ant_norm += segreg_calculator::mlm_ant_vec[i];
     pos_norm += segreg_calculator::mlm_pos_vec[i];
 }

 for (unsigned i = 0; i != segreg_calculator::mlm_ant_vec.size(); i++) { // now rescale
    segreg_calculator::mlm_ant_vec[i] = segreg_calculator::mlm_ant_vec[i]/ant_norm;
     segreg_calculator::mlm_pos_vec[i] = segreg_calculator::mlm_pos_vec[i]/pos_norm;
 }

for (unsigned i = 0; i != segreg_calculator::mlm_ant_vec.size(); i++) { // the likelihood after rescaling
    rescale_like += segreg_calculator::mlm_ant_vec[i]*segreg_calculator::mlm_pos_vec[i]; 
 }


 log_double new_log_like(0.0);
 new_log_like = ant_norm*pos_norm*rescale_like; // due to JA
 
  delete mpl;

//  cout << ps->name() << ' ' << like << endl;


//  return like;
 new_log_like = ant_norm*pos_norm*rescale_like; // due to JA
    return new_log_like;
}


//----------------------------------------------------------------------------
//
//      calc_unconnected_regressive(...)
//
//----------------------------------------------------------------------------           

void segreg_calculator::calc_unconnected_regressive()
{
  for(PedigreeDataSet::MemberIterator
          unconnected = my_ped_data.get_unconnected_begin();
          unconnected != my_ped_data.get_unconnected_end(); ++unconnected)
  {
    log_double like = calc_unconnected_regressive(*unconnected, false);

    if(my_asc_like_elts)
      like /= calc_unconnected_regressive(*unconnected, true);

    if(!accumulate_likelihood(like)) return;
  }
}

//----------------------------------------------------------------------------
//
//      calc_unconnected_FPMM(...)
//
//----------------------------------------------------------------------------           

void segreg_calculator::calc_unconnected_FPMM()
{
  for(PedigreeDataSet::MemberIterator
          unconnected = my_ped_data.get_unconnected_begin();
          unconnected != my_ped_data.get_unconnected_end(); ++unconnected)
  {
    log_double like = calc_unconnected_FPMM(*unconnected, false);

    if(asc_fpmmsl)
      like /= calc_unconnected_FPMM(*unconnected, true);

    if(!accumulate_likelihood(like)) return;
  }
}

//----------------------------------------------------------------------------
//
//      calc_unconnected_MLM(...)
//
//----------------------------------------------------------------------------           

void segreg_calculator::calc_unconnected_MLM()
{
  for(PedigreeDataSet::MemberIterator
          unconnected = my_ped_data.get_unconnected_begin();
          unconnected != my_ped_data.get_unconnected_end(); ++unconnected)
  {
    log_double like = calc_unconnected_MLM(*unconnected, false);

    if(abmc)
      like /= calc_unconnected_MLM(*unconnected, true);

    if(!accumulate_likelihood(like)) return;
  }
}

//----------------------------------------------------------------------------------
//
// calc_unconnected_regressive(...)
//
//----------------------------------------------------------------------------------
/// unconnected likelihood calculation for a single member

inline
log_double
segreg_calculator::calc_unconnected_regressive
   (const FPED::Member& member, bool asc)
{
  log_double like(0.0);

  RegUnconnectedLikelihood p(((asc) ? my_asc_like_elts->get_pc() : my_like_elts->get_pc()), md);
    
  for(int genotype = 0; genotype < 3; ++genotype)
  { 
    like += p.unconnected_likelihood(genotype, &member);
  }
  
  return like;
}      

//----------------------------------------------------------------------------------
//
// calc_unconnected_FPMM(...)
//
//----------------------------------------------------------------------------------
/// unconnected likelihood calculation for a single member

inline
log_double
segreg_calculator::calc_unconnected_FPMM
   (const FPED::Member& member, bool asc)
{
  log_double like(0.0);

  FPMM_SL* p = (asc) ? asc_fpmmsl : fpmmsl;
    
  for(genotype_index geno = index_AA; geno != index_INVALID; ++geno)
  {
    for(size_t poly = 0; poly < md.fpmm_sub_model.max_pgt(); ++poly)
    { 
      like += p->unconnected_likelihood(&member, geno, poly);
    }
  }
  
  return like;
}      

//----------------------------------------------------------------------------------
//
// calc_unconnected_MLM(...)
//
//----------------------------------------------------------------------------------
/// unconnected likelihood calculation for a single member

inline
log_double
segreg_calculator::calc_unconnected_MLM
   (const FPED::Member& member, bool asc)
{
  log_double like(0.0);

  binary_member_calculator* p = (asc) ? abmc : bmc;
    
  for(genotype_index genotype = index_AA; genotype != index_INVALID; ++genotype)
  { 
    like += p->unconnected_likelihood(member, genotype);
  }
  
  return like;
}      

void
segreg_calculator::calc_connected_pgeno_regressive
    (const FPED::Subpedigree& s, pg::post_geno_map& pmap)
{
  // Create the peeler

  regressive_peeler rpl(s, *my_like_elts);

  for(size_t i = 0; i < s.member_count(); ++i)
  {
    calc_member_pgeno_regressive(pmap, rpl, s.member_index(i));
  }
}

void
segreg_calculator::calc_connected_pen_func
    (const FPED::Subpedigree& s, pf::pen_func_map& pmap)
{
  // Create the peeler

  regressive_peeler rpl(s, *my_like_elts);

  for(size_t i = 0; i < s.member_count(); ++i)
  {
    if(MPED::mp_utilities::is_founder(s.member_index(i)))
      calc_founder_pen_func(s.member_index(i), pmap);
    else
      calc_nonfounder_pen_func(s.member_index(i), pmap);
  }
}

void
segreg_calculator::calc_connected_pgeno_FPMM
    (const FPED::Subpedigree& s, pg::post_geno_map& pmap)
{
  // Create the peeler

  FPMM_peeler fpl(s, md);

  fpl.set_SL(fpmmsl);

  fpmmsl->set_peeler(&fpl);

  for(size_t i = 0; i < s.member_count(); ++i)
  {
    calc_member_pgeno_FPMM(pmap, s.member_index(i));
  }
}

void
segreg_calculator::calc_connected_pgeno_MLM
    (const FPED::Subpedigree& s, pg::post_geno_map& pmap)
{
  // Create the peeler

  MlmLikelihoodElements lelt(*my_ped_data.get_raw_data(),md,false);
  mpl = new mlm_peeler(s,bmc,*my_ped_data.get_raw_data(),md,lelt,false);
                       

  for(size_t i = 0; i < s.member_count(); ++i)
  {
    calc_member_pgeno_MLM(pmap, s.member_index(i));
  }

  delete mpl;
}

void segreg_calculator::calc_unconnected_pgeno_regressive
  (const FPED::Member& m, pg::post_geno_map& pmap)
{
  // Calculate the posterior probabilities of the genotypes

  RegUnconnectedLikelihood p(my_like_elts->get_pc(), md);

  log_double like_AA = p.unconnected_likelihood(index_AA, &m);
  log_double like_AB = p.unconnected_likelihood(index_AB, &m);
  log_double like_BB = p.unconnected_likelihood(index_BB, &m);

  double prob_AA, prob_AB, prob_BB;

  renormalize(prob_AA, prob_AB, prob_BB, like_AA, like_AB, like_BB);

  // Get the penetrances
  const TypeDescription& tdesc = my_like_elts->get_type_description();

  double pen_AA = my_like_elts->get_penetrance(m, tdesc.get_state(index_AA));
  double pen_AB = my_like_elts->get_penetrance(m, tdesc.get_state(index_AB));
  double pen_BB = my_like_elts->get_penetrance(m, tdesc.get_state(index_BB));

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));
}

void segreg_calculator::calc_unconnected_pen_func
  (const FPED::Member& m, pf::pen_func_map& pmap)
{
  // Get the penetrances
  const TypeDescription& tdesc = my_like_elts->get_type_description();

  double pen_AA = my_like_elts->get_penetrance(m, tdesc.get_state(index_AA));
  double pen_AB = my_like_elts->get_penetrance(m, tdesc.get_state(index_AB));
  double pen_BB = my_like_elts->get_penetrance(m, tdesc.get_state(index_BB));

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(pen_AA) pmap[&m].push_back(pf::pen_type("A/A", pen_AA));
  if(pen_AB) pmap[&m].push_back(pf::pen_type("A/B", pen_AB));
  if(pen_BB) pmap[&m].push_back(pf::pen_type("B/B", pen_BB));
}

void segreg_calculator::calc_unconnected_pgeno_FPMM
(const FPED::Member&  m, pg::post_geno_map&  pmap)
{
  log_double like_AA(0.0);
  log_double like_AB(0.0);
  log_double like_BB(0.0);

  double prob_AA, prob_AB, prob_BB;

  double pen_AA = 0.0;
  double pen_AB = 0.0;
  double pen_BB = 0.0;

  for(size_t i = 0; i <= md.fpmm_sub_model.loci() * 2; ++i)
  {
    like_AA += fpmmsl->unconnected_likelihood(&m, index_AA, i);
    like_AB += fpmmsl->unconnected_likelihood(&m, index_AB, i);
    like_BB += fpmmsl->unconnected_likelihood(&m, index_BB, i);

    pen_AA += fpmmsl->penetrance(&m, index_AA, i);
    pen_AB += fpmmsl->penetrance(&m, index_AB, i);
    pen_BB += fpmmsl->penetrance(&m, index_BB, i);
  }

  renormalize(prob_AA, prob_AB, prob_BB, like_AA, like_AB, like_BB);

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));
}

void segreg_calculator::calc_unconnected_pgeno_MLM
  (const FPED::Member&  m, pg::post_geno_map&  pmap)
{
  double prob_AA = bmc->unconnected_likelihood(m, index_AA);
  double prob_AB = bmc->unconnected_likelihood(m, index_AB);
  double prob_BB = bmc->unconnected_likelihood(m, index_BB);

  renormalize(prob_AA, prob_AB, prob_BB);

  // Get the penetrances

  double pen_AA = bmc->get_penetrance(m, index_AA);
  double pen_AB = bmc->get_penetrance(m, index_AB);
  double pen_BB = bmc->get_penetrance(m, index_BB);

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));
}

void segreg_calculator::calc_member_pgeno_regressive
(pg::post_geno_map& pmap,
 regressive_peeler& rp,
 const FPED::Member& m)
{
  log_double like_AA, like_AB, like_BB;

  // Get the posterior genotype probabilities
  const TypeDescription& tdesc = my_like_elts->get_type_description();

  like_AA = calc_subped_given_member_regressive(&rp, m, tdesc.get_state(index_AA));
  like_AB = calc_subped_given_member_regressive(&rp, m, tdesc.get_state(index_AB));
  like_BB = calc_subped_given_member_regressive(&rp, m, tdesc.get_state(index_BB));

  double prob_AA, prob_AB, prob_BB;

  renormalize(prob_AA, prob_AB, prob_BB, like_AA, like_AB, like_BB);

  // Get the penetrances

  double pen_AA = my_like_elts->get_penetrance(m, tdesc.get_state(index_AA));
  double pen_AB = my_like_elts->get_penetrance(m, tdesc.get_state(index_AB));
  double pen_BB = my_like_elts->get_penetrance(m, tdesc.get_state(index_BB));

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store all values > cutoff

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));
}

void segreg_calculator::calc_member_pgeno_MLM
(pg::post_geno_map& pmap, const FPED::Member& m)
{
  log_double like_AA, like_AB, like_BB;

  MlmLikelihoodElements lelt(*my_ped_data.get_raw_data(),md,false);
  
  TypeDescription::StateIterator state = lelt.get_type_description().begin();

  like_AA = calc_subped_given_member_MLM(mpl, m, *state); ++state;
  like_AB = calc_subped_given_member_MLM(mpl, m, *state); ++state;
  like_BB = calc_subped_given_member_MLM(mpl, m, *state); ++state;

  double prob_AA, prob_AB, prob_BB;

  renormalize(prob_AA, prob_AB, prob_BB, like_AA, like_AB, like_BB);
  
  // Get the penetrances

  double pen_AA = bmc->get_penetrance(m, index_AA);
  double pen_AB = bmc->get_penetrance(m, index_AB);
  double pen_BB = bmc->get_penetrance(m, index_BB);

  renormalize(pen_AA, pen_AB, pen_BB);
  
  // Store all values > cutoff

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));

}

void segreg_calculator::calc_member_pgeno_FPMM
(pg::post_geno_map& pmap, const FPED::Member& m)
{
  log_double like_AA(0.0);
  log_double like_AB(0.0);
  log_double like_BB(0.0);

  double prob_AA, prob_AB, prob_BB;

  double pen_AA = 0.0;
  double pen_AB = 0.0;
  double pen_BB = 0.0;

  FPMM_peeler* pl = (FPMM_peeler*) fpmmsl->get_peeler();

  for(size_t i = 0; i <= md.fpmm_sub_model.loci() * 2; ++i)
  {
    like_AA += calc_subped_given_member_FPMM(pl, m, genetic_info(index_AA, i));
    like_AB += calc_subped_given_member_FPMM(pl, m, genetic_info(index_AB, i));
    like_BB += calc_subped_given_member_FPMM(pl, m, genetic_info(index_BB, i));

    pen_AA += fpmmsl->penetrance(&m, index_AA, i);
    pen_AB += fpmmsl->penetrance(&m, index_AB, i);
    pen_BB += fpmmsl->penetrance(&m, index_BB, i);
  }

  renormalize(prob_AA, prob_AB, prob_BB, like_AA, like_AB, like_BB);

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(prob_AA || pen_AA) pmap[&m].push_back(pg::pen_type("A/A", prob_AA, pen_AA));
  if(prob_AB || pen_AB) pmap[&m].push_back(pg::pen_type("A/B", prob_AB, pen_AB));
  if(prob_BB || pen_BB) pmap[&m].push_back(pg::pen_type("B/B", prob_BB, pen_BB));
}

void
segreg_calculator::calc_nonfounder_pen_func
    (const FPED::Member& ind, pf::pen_func_map& pmap)
{
  const TypeDescription& tdesc = my_like_elts->get_type_description();

  PenetranceContext context = my_like_elts->get_penetrance_context();

  // Add the nuclear family to the context, but exclude siblings.
  context.set_nuclear_family(*ind.family(), false);

  // Create a vector to store the child's likelihoods in.  This should be 27
  // (genotypes^3) elements per child, of which many will be zero, but for
  // simplicity, the memory waste is acceptable.

  vector<log_double> child_likelihoods(27, log_double(0.0));

  for(TypeDescription::StateIterator mgeno = tdesc.begin();
      mgeno != tdesc.end(); ++mgeno)
  {
    context.set_mother_state(*mgeno);
    for(TypeDescription::StateIterator fgeno = tdesc.begin();
        fgeno != tdesc.end(); ++fgeno)
    {
      context.set_father_state(*fgeno);

      double trans[3];
      
      trans[index_AA] = md.transm_sub_model.prob(index_AA, mgeno->get_index(), fgeno->get_index());
      trans[index_AB] = md.transm_sub_model.prob(index_AB, mgeno->get_index(), fgeno->get_index());
      trans[index_BB] = md.transm_sub_model.prob(index_BB, mgeno->get_index(), fgeno->get_index());
      
      // The index into the 27 element matrix is 9 * mother + 3 * father + child

      size_t mfindex = mgeno->get_index() * 9 + fgeno->get_index() *3;

      // For each child state

      for(TypeDescription::StateIterator cgeno = tdesc.begin();
          cgeno != tdesc.end(); ++cgeno)
      {

        // If the transmission is zero, this triple is not possible, and we
        // need do no work for it.
        if(trans[cgeno->get_index()] == 0.0) continue;

        child_likelihoods[mfindex + cgeno->get_index()] =
            my_like_elts->get_penetrance(ind, *cgeno, context);
      }
    }
  }
  
  // Normalize the child to 1 and start inserting them into the
  // penetrance_map

  log_double child_sum(0.0);

  for(size_t i = 0; i < 27; ++i)
  {
    child_sum += child_likelihoods[i];
  }

  for(size_t i = 0; i < 27; ++i)
  {
    string geno = "";
    string mgeno = "";
    string fgeno = "";
    
    log_double d = child_likelihoods[i] / child_sum;
    
    if(d.get_double() < PEN_CUTOFF) continue;

    // Mother state
    
    if     (i <  9)  mgeno ="A/A";
    else if(i < 18)  mgeno ="A/B";
    else             mgeno ="B/B";
    
    // Father state
    
    if     (i%9 < 3)  fgeno ="A/A";
    else if(i%9 < 6)  fgeno ="A/B";
    else              fgeno ="B/B";

    // Child state
    
    if     (i%3 == 0)  geno ="A/A";
    else if(i%3 == 1)  geno ="A/B";
    else               geno ="B/B";
    
    pmap[&ind].push_back(pf::pen_type(geno, mgeno, fgeno, d.get_double()));
  }
}

void segreg_calculator::calc_founder_pen_func
  (const FPED::Member& m, pf::pen_func_map& pmap)
{
  const TypeDescription& tdesc = my_like_elts->get_type_description();

  // Get the penetrances
  double pen_AA = my_like_elts->get_penetrance(m, tdesc.get_state(index_AA));
  double pen_AB = my_like_elts->get_penetrance(m, tdesc.get_state(index_AB));
  double pen_BB = my_like_elts->get_penetrance(m, tdesc.get_state(index_BB));

  renormalize(pen_AA, pen_AB, pen_BB);

  // Store non zero values.

  if(pen_AA) pmap[&m].push_back(pf::pen_type("A/A", pen_AA));
  if(pen_AB) pmap[&m].push_back(pf::pen_type("A/B", pen_AB));
  if(pen_BB) pmap[&m].push_back(pf::pen_type("B/B", pen_BB));
}

bool
segreg_calculator::accumulate_likelihood
    (log_double likelihood)
{
    likelihood1 *= likelihood;

    return finite(likelihood1.get_log());
}

bool segreg_calculator::finalize_likelihood()
{
  return finite(likelihood1.get_log());
}

}}
