#include "segreg/freq_estimates.h"

namespace SAGE
{
namespace SEGREG
{

//
// option refers to hwe, nhwe wtc.
// in the set calls the q_A corresponds to the frequence of A
// the next three arguments refer to genotype probabilities
// 

bool genotype_frequency_initial_estimates::cont;
 
bool freq_model_valid(const genotype_frequency_sub_model& m)
{
  return !SAGE::isnan(m.freq_A())       &&
         !SAGE::isnan(m.prob(index_AA)) &&
         !SAGE::isnan(m.prob(index_AB)) &&
         !SAGE::isnan(m.prob(index_BB));
}
/* +++++ Original version here 
genotype_frequency_initial_estimates::genotype_frequency_initial_estimates
    (const mean_split& splitter, const freq_sub_model& mod)
{
  // Get the model option

  genotype_frequency_sub_model::sm_option option = mod.option();

  //lint --e{534}  We ignore lots of returns here

  // Check the model for correctness.  If correct, just copy it.

  if(option != genotype_frequency_sub_model::NONE && freq_model_valid(mod))
  {
    my_model_count = 3;

    my_estimates[0] = mod;
    my_estimates[1] = mod;
    my_estimates[2] = mod;

    return;
  }

  // If the option was NONE, we assume hwe

  if(option == genotype_frequency_sub_model::NONE)
    option = genotype_frequency_sub_model::hwe;

  // Generate our estimates.  These numbers are based upon the SEGREG
  // formula documentation, Appendix B.

  // Get the number of non-missing values from the mean_split object

  size_t trait_count = splitter.get_n();

  // Set first estimate

  double qA = 1.0 - sqrt(1.0 - 1.0/trait_count);

  my_estimates[0].set
      (option, qA, QNAN, QNAN, QNAN, model_input(0.0, true), false);

  // Set second estimate

  qA = sqrt(1.0 - 1.0/trait_count);

  my_estimates[1].set
      (option, qA, QNAN, QNAN, QNAN, model_input(0.0, true), false);

  // Set third estimate if there is a valid quantity for it.

  bool b = splitter.get_status();

  if(!b)
  {
    my_model_count = 2;
  }
  else
  {
    my_model_count = 3;

    
    qA = 1.0 - sqrt(1.0 - splitter.get_p());

    my_estimates[2].set
        (option, qA, QNAN, QNAN, QNAN,  model_input(0.0, true), false);
  }
}
*/
//
// A modified method which takes into account the discrepancy in 
// frequency estimates between the discrete and continuous cases
//

genotype_frequency_initial_estimates::genotype_frequency_initial_estimates
    (const mean_split& splitter, const freq_sub_model& mod)
{
  // Get the model option

  genotype_frequency_sub_model::sm_option option = mod.option();

  //lint --e{534}  We ignore lots of returns here

  // Check the model for correctness.  If correct, just copy it.

  if(option != genotype_frequency_sub_model::NONE && freq_model_valid(mod))
  {
    my_model_count = 3;

    my_estimates[0] = mod;
    my_estimates[1] = mod;
    my_estimates[2] = mod;

   return;
  }

  // If the option was NONE, we assume hwe

  if(option == genotype_frequency_sub_model::NONE)
    option = genotype_frequency_sub_model::hwe;

  // Get the number of non-missing values from the mean_split object

  size_t trait_count = splitter.get_n();

  // Set first estimate

  if (genotype_frequency_initial_estimates::cont){

 double qA = 1.0 - sqrt(1.0 - 1.0/trait_count);

  my_estimates[0].set
      (option, qA, QNAN, QNAN, QNAN, model_input(0.0, true), false);

  // Set second estimate

  qA = sqrt(1.0 - 1.0/trait_count);

  my_estimates[1].set
      (option, qA, QNAN, QNAN, QNAN, model_input(0.0, true), false);

  // Set third estimate if there is a valid quantity for it.

  bool b = splitter.get_status();

  if(!b)
  {
    my_model_count = 2;
  }
  else
  {
    my_model_count = 3;

    
    qA = 1.0 - sqrt(1.0 - splitter.get_p());

    my_estimates[2].set
        (option, qA, QNAN, QNAN, QNAN,  model_input(0.0, true), false);
  }
 } // trait continuous 
 else{ // binary trait
       
    for(unsigned i = 0; i != 5; ++i){
     double p;
       if (i == 0) p = 1.0/(trait_count);    
       if ((i > 0) && ( i < 5)) p = (i*1.0)/5.0;
       if (i == 5) p = (trait_count -1.0)/trait_count;
       double qA = 1.0 - sqrt(1.0 - p);
       my_estimates[i].set\
        (option, qA, QNAN, QNAN, QNAN,  model_input(0.0, true), false);
     }
    my_model_count = 6;
  } 
}

}
}
