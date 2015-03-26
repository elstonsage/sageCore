
//==================================================================================
// File:        segreg_calculator.ipp
//
// Purpose:     inline functions for construction and implementing calculation for
//              segregation likelihood.
//              When constructing object, it checks the validation of the pedigree.
//
// Author:      Kai He
//
// History:     07/26/2001 initialed
//
// Copyright (c) 2001 R. C. Elston
//==================================================================================


#ifndef SEGREG_CALCULATOR_H
#include "segreg/segreg_calculator.h"
#endif

namespace SAGE
{
namespace SEGREG
{

//----------------------------------------------------------------------------------
//
// Construction
//
//----------------------------------------------------------------------------------

inline
segreg_calculator::segreg_calculator
    (const PedigreeDataSet& ped_data,
     const model&           m)
  : my_ped_data    (ped_data),
    md             ( m  ),

    my_like_elts   (NULL),
    my_asc_like_elts(NULL),
    fpmmsl         (NULL),
    asc_fpmmsl     (NULL),
    bmc            (NULL),
    abmc           (NULL),
    my_mlm_corr_verifier (NULL),
    last_likelihood(QNAN)
{ 
   nfe = 0;

   setup_components(); 

   set_continuous_penalty_component (1.0   );
}

//----------------------------------------------------------------------------------
//
// Destruction 
//
//----------------------------------------------------------------------------------

inline segreg_calculator::~segreg_calculator()
{
  if(my_like_elts) { delete my_like_elts; my_like_elts = NULL; }

  if(fpmmsl)     { delete fpmmsl;     fpmmsl     = NULL; }

  if(my_asc_like_elts)    { delete my_asc_like_elts; my_asc_like_elts = NULL; }
  if(asc_fpmmsl) { delete asc_fpmmsl; asc_fpmmsl = NULL; }

  if(bmc)        { delete bmc;        bmc        = NULL; }
  if(abmc)       { delete abmc;       abmc       = NULL; }

  my_mlm_corr_verifier = std::auto_ptr<MlmCorrelationVerifier>(NULL);

  //if(mlm_resid_corr) { delete mlm_resid_corr; mlm_resid_corr = NULL; }
}

//----------------------------------------------------------------------------------
//
// setup_components(...)
//
//----------------------------------------------------------------------------------

inline
void segreg_calculator::setup_components()
{
   model_class mc = md.get_model_class();
 
   switch(mc)
   {
     case model_A    : 
     case model_D    : setup_regressive_components(); return;
     case model_FPMM : setup_fpmm_components(); return;
     case model_MLM  : setup_mlm_components(); return;

     default         : break;
   }
}

inline
void segreg_calculator::setup_regressive_components()
{
  my_like_elts = new LikelihoodElements(*my_ped_data.get_raw_data(), md, false);
  if(using_ascertainment())
  {
    my_asc_like_elts = new LikelihoodElements(*my_ped_data.get_raw_data(), md, true);
  }
}

inline
void segreg_calculator::setup_mlm_components()
{
  bmc = new binary_member_calculator(*my_ped_data.get_raw_data(),md,false);
  
  my_mlm_corr_verifier =
    std::auto_ptr<MlmCorrelationVerifier>
        (new MlmCorrelationVerifier(my_ped_data.get_subpedigrees(),
                                    md.resid_sub_model, *bmc));
  
  if(using_ascertainment())
  {
    abmc = new binary_member_calculator(*my_ped_data.get_raw_data(),md,true);
  }
}

inline
void segreg_calculator::setup_fpmm_components()
{
  fpmmsl = new FPMM_SL(*my_ped_data.get_raw_data(), md, false);

  if(using_ascertainment())
  {
    asc_fpmmsl = new FPMM_SL(*my_ped_data.get_raw_data(), md, true);            
  }
}

inline
bool
  segreg_calculator::using_ascertainment() const
{
  return md.ascer_sub_model.s_option() != ascertainment_sub_model::none;
}

//----------------------------------------------------------------------------------
//
// set penalization component - C 
//
//----------------------------------------------------------------------------------
    
inline       
void segreg_calculator::set_continuous_penalty_component(double c)
{
   c_penalty_value = max(c, 0.0);
}

inline
double segreg_calculator::calculate_continuous_penalty() const
{
   double continuous_penalty = 0.0;

   // Get Psi values

   double PSI_AA, PSI_AB, PSI_BB;

   PSI_AA  = md.freq_sub_model.prob(index_AA);
   PSI_AB  = md.freq_sub_model.prob(index_AB);          
   PSI_BB  = md.freq_sub_model.prob(index_BB);          

   // Calculate based on mean model

   typedef genotype_specific_mean_sub_model	msm;

   msm::sm_option sp = md.mean_sub_model.option();

   double c = c_penalty_value;

   switch(sp)
   {
     case msm::one       : break;
     case msm::two       : 
     case msm::two_dom   : continuous_penalty = c * log(4.0 * PSI_BB * (1.0 - PSI_BB)); break;
     case msm::two_rec   : continuous_penalty = c * log(4.0 * PSI_AA * (1.0 - PSI_AA)); break;
     case msm::three     : 
     case msm::three_add : 
     case msm::three_dec : 
     case msm::three_inc : continuous_penalty = c * log(8.0 * PSI_AA * PSI_AB * PSI_BB); break;
     default             : break;
   }

   return continuous_penalty;
}    
 
//----------------------------------------------------------------------------------
//
// set_discrete_penalty_component(...)
//
// calculate population prevalence: Prev
//----------------------------------------------------------------------------------

inline
log_double segreg_calculator::calculate_prevalence_penalty() const
{
  return md.prev_sub_model.get_prevalence_penalty();
}

//----------------------------------------------------------------------------------
//
// get_subpedigree_count(...)
//
//----------------------------------------------------------------------------------

inline
uint segreg_calculator::get_subpedigree_count() const
{
   return my_ped_data.get_subpedigree_count();
} 
//----------------------------------------------------------------------------------
//
// get_unconnected_count(...)
//
//----------------------------------------------------------------------------------
  
inline  
uint segreg_calculator::get_unconnected_count() const
{   
   return my_ped_data.get_unconnected_count();                 
}   

/// Calculates the total likelihood by calculating the connected and
/// unconnected portions and multiplying

inline
log_double segreg_calculator::calc_likelihood()
{
  clear_likelihood();

  model_class  mc = md.get_model_class();

  switch(mc)
  {
    case model_A    :
    case model_D    :
      calc_connected_regressive   ();
      calc_unconnected_regressive ();
      break;

    case model_MLM  :
      calc_connected_MLM          ();
      calc_unconnected_MLM        ();
      break;

    case model_FPMM :
      calc_connected_FPMM         ();
      calc_unconnected_FPMM       ();
      break;

    default         :
      SAGE_internal_error();
      break;
  }

  finalize_likelihood();

  return likelihood1;
}

//----------------------------------------------------------------------------------
//
// calc_subped_given_member_regressive(...) 
//
//----------------------------------------------------------------------------------
// Subpedigree probability calculation for Regressive Model
// This function will go outside to wake up  regressive_peeler

inline
log_double segreg_calculator::calc_subped_given_member_regressive
(regressive_peeler* rpl, const member_type& indi, const TypeDescription::State&  genotype)
{
   log_double reg_prob(0.0);
   
   log_double ant_prob = rpl->anterior(indi, genotype);
   
   log_double pos_prob = rpl->posterior(indi, genotype);

   segreg_calculator::reg_ant_vec.push_back(ant_prob.get_double()); // due to JA
   segreg_calculator::reg_pos_vec.push_back(pos_prob.get_double()); // due to JA

   reg_prob = ant_prob * pos_prob;

   return reg_prob;
}

inline
log_double segreg_calculator::calc_subped_given_member_MLM
(mlm_peeler*  mpl, const member_type& i, const TypeDescription::State& g)
{
//   cout << "Inside ant-pos part " << endl; // for debugging purposes (due to JA)
   log_double mlm_prob(0.0);

   log_double ant_prob = mpl->anterior (i, g);
   log_double pos_prob = mpl->posterior(i, g);

   segreg_calculator::mlm_ant_vec.push_back(ant_prob.get_double()); // due to JA
   segreg_calculator::mlm_pos_vec.push_back(pos_prob.get_double()); // due to JA

   mlm_prob = ant_prob * pos_prob;

   return mlm_prob;
}

//----------------------------------------------------------------------------------
// 
// calc_subped_given_member_FPMM(...)
//
//----------------------------------------------------------------------------------

inline 
log_double segreg_calculator::calc_subped_given_member_FPMM
(FPMM_peeler* fpl, const member_type& indi, genetic_info gi)
{
   log_double fpmm_prob (0.0);
   
   log_double ant_prob = fpl->anterior (indi, gi);
  
   log_double pos_prob = fpl->posterior(indi, gi);
  
   segreg_calculator::fpmm_ant_vec.push_back(ant_prob.get_double()); // due to JA
   segreg_calculator::fpmm_pos_vec.push_back(pos_prob.get_double()); // due to JA

   fpmm_prob = ant_prob * pos_prob;
  
   return fpmm_prob;

}

inline void segreg_calculator::calculate_pen_func_probs(pf::pen_func_map& p)
{
  model_class  mc = md.get_model_class();

  if(mc != model_A && mc != model_D) return;

  for(PedigreeDataSet::SubpedigreeIterator
          ps_itr = my_ped_data.get_subpedigree_begin();
          ps_itr != my_ped_data.get_subpedigree_end(); ++ps_itr)
  {
    calc_connected_pen_func(*ps_itr, p);
  }

  for(PedigreeDataSet::MemberIterator
          mem_itr = my_ped_data.get_unconnected_begin();
          mem_itr != my_ped_data.get_unconnected_end(); ++mem_itr)
  {
    calc_unconnected_pen_func(*mem_itr, p);
  }

}

inline void segreg_calculator::calculate_post_geno_probs(pg::post_geno_map& p)
{
  void (segreg_calculator::* sped_func) (const FPED::Subpedigree&,  pg::post_geno_map&);
  void (segreg_calculator::* mem_func)  (const FPED::Member&,       pg::post_geno_map&);

  model_class  mc = md.get_model_class();

  switch(mc)
  {
    case model_A    :
    case model_D    : sped_func = &segreg_calculator::calc_connected_pgeno_regressive;
                      mem_func  = &segreg_calculator::calc_unconnected_pgeno_regressive;
                      break;
    case model_MLM  : sped_func = &segreg_calculator::calc_connected_pgeno_MLM;
                      mem_func  = &segreg_calculator::calc_unconnected_pgeno_MLM;
                      break;
    case model_FPMM : sped_func = &segreg_calculator::calc_connected_pgeno_FPMM;
                      mem_func  = &segreg_calculator::calc_unconnected_pgeno_FPMM;
                      break;
    case model_INVALID : // Will never happen
                      break;
  }

  for(PedigreeDataSet::SubpedigreeIterator
          ps_itr = my_ped_data.get_subpedigree_begin();
          ps_itr != my_ped_data.get_subpedigree_end(); ++ps_itr)
  {
    (this->*sped_func)(*ps_itr, p);
  }

  for(PedigreeDataSet::MemberIterator
          unconnected = my_ped_data.get_unconnected_begin();
          unconnected != my_ped_data.get_unconnected_end(); ++unconnected)
  {
    (this->*mem_func)(*unconnected, p);
  }
}

inline void
segreg_calculator::calculate_mlm_resid_corr
    (MlmResidCorrelationCalculator&  mra) const
{
   mra.calculate(*bmc, my_ped_data.get_subpedigrees(),md);
}

inline void segreg_calculator::clear_likelihood()
{
  likelihood1 = 1.0;
  likelihood2 = 1.0;
}

//---------------------------------------------------------------------------------

}}
