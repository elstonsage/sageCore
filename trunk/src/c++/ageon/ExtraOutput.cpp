#include <string>
#include <fstream>
#include "numerics/functions.h"
#include "ageon/ExtraOutput.h"

namespace SAGE {
namespace AO   {

//===================================
//
//  generateFile(...)
//
//===================================
void
ExtraOutput::generateFile(
  const SAMPLING::PartitionedMemberDataSample & sample,
  const ModelTraitsVector                     & model_traits, 
  const std::string                           & analysis_name,
        std::ostream                          & info)
{
  // Grab the individual / sex codes:
  const std::string & ind_missing = sample.getMultipedigree().info().individual_missing_code (),
                    & sex_missing = sample.getMultipedigree().info().sex_code_unknown        (),
                    & sex_male    = sample.getMultipedigree().info().sex_code_male           (),
                    & sex_female  = sample.getMultipedigree().info().sex_code_female         (),
                      filename    = analysis_name;

  // Pedigree block info for .inf file:
  info << endl
       << "Note:\n"
       << "  To use the susceptability traits and survival analysis residuals,\n"
       << "  use the additinal output files from this analysis\n"
       << "  '" << analysis_name+".par' and '" << analysis_name+".ped'.\n\n";

  // Parameter file:
  std::ofstream opar((filename+".par").c_str());

  opar << "pedigree, file=\"" << (filename+".ped") << "\"\n"
       << "{\n"
       << "  delimiter_mode = multiple\n"
       << "  delimiters=\" \t\"\n"
       << "  individual_missing_value=\"" << ind_missing << "\"\n"
       << "  sex_code, male=\"" << sex_male << "\", female=\"" << sex_female << "\", unknown=\"" << sex_missing << "\"\n"
       << "  pedigree_id   = PED\n"
       << "  individual_id = ID\n"
       << "  parent_id     = P1\n"
       << "  parent_id     = P2\n"
       << "  sex_field     = SEX\n"
       << "  trait = nt_equal_trait, missing=\"nan\"\n"
       << "  trait = nt_equal_residual, missing=\"nan\"\n"
       << "  trait = nt_free_trait, missing=\"nan\"\n"
       << "  trait = nt_free_residual, missing=\"nan\"\n"
       << "  trait = t_equal_trait, missing=\"nan\"\n"
       << "  trait = t_equal_residual, missing=\"nan\"\n"
       << "  trait = t_free_trait, missing=\"nan\"\n"
       << "  trait = t_free_residual, missing=\"nan\"\n"
       << "}\n\n"
       << std::endl
       << std::flush;

  // Pedigree file:
  std::ofstream o((filename+".ped").c_str());

  o << "PED\tID\tP1\tP2\tSEX\t"
    << "nt_equal_trait\t"
    << "nt_equal_residual\t"
    << "nt_free_trait\t"
    << "nt_free_residual\t"
    << "t_equal_trait\t"
    << "t_equal_residual\t"
    << "t_free_trait\t"
    << "t_free_residual\n";
    
  for(size_t i = 0; i < sample.getMultipedigree().member_count(); ++i)
  {
    // Grab the member reference:
    const FPED::Member & member = sample.getMultipedigree().member_index(i);

    // Grab sex:
    bool male   = member.is_male();
    bool female = member.is_female();
    
    // Structure:
    o << member.pedigree()->name()                                   << "\t"
      << member.name()                                               << "\t"
      << (member.parent1() ? member.parent1()->name() : ind_missing) << "\t"
      << (member.parent2() ? member.parent2()->name() : ind_missing) << "\t";
    
    // Sex:
    o << (male ? sex_male : (female ? sex_female : sex_missing)) << "\t";
    
    // Traits:
    o << model_traits[NO_TRUNCATION  | SUSCEPTIBILITIES_EQUAL] . cond_probs[i] << "\t"
      << model_traits[NO_TRUNCATION  | SUSCEPTIBILITIES_EQUAL] . surv_resid[i] << "\t"
      << model_traits[NO_TRUNCATION  | SUSCEPTIBILITIES_FREE]  . cond_probs[i] << "\t"
      << model_traits[NO_TRUNCATION  | SUSCEPTIBILITIES_FREE]  . surv_resid[i] << "\t"
      << model_traits[USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL] . cond_probs[i] << "\t"
      << model_traits[USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL] . surv_resid[i] << "\t"
      << model_traits[USE_TRUNCATION | SUSCEPTIBILITIES_FREE]  . cond_probs[i] << "\t"
      << model_traits[USE_TRUNCATION | SUSCEPTIBILITIES_FREE]  . surv_resid[i];
      
    // Newline:
    o << std::endl;
  }
}

//===================================
//
//  populateModelTraits(...)
//
//===================================
void 
ExtraOutput::populateModelTraits(const MemberCovariateCalculator & mcc, ModelTraits & model_traits)
{
  model_traits.cond_probs.resize(mcc.getSample().getMultipedigree().member_count(), 0.0);
  model_traits.surv_resid.resize(mcc.getSample().getMultipedigree().member_count(), 0.0);
  
  for(size_t i = 0; i < mcc.getSample().getMultipedigree().member_count(); ++i)
  {
    bool   affected       = mcc.getSample().getAdjValue(i, "CORE_TRAITS", "affectedness") == 1;
    double suscept        = mcc.get_genetic_suscept (i),
           AO_transf      = mcc.get_AO_transf       (i),
           AE_transf      = mcc.get_AE_transf       (i),
           AO_mean        = mcc.get_AO_mean         (i),
           AO_stdev       = mcc.get_AO_stdev        (i),
           aff_cdf_term   = normal_cdf((AO_transf - AO_mean) / AO_stdev),
           unaff_cdf_term = normal_cdf((AE_transf - AO_mean) / AO_stdev);
           
    if(affected)
    {
      model_traits.cond_probs[i] = 1.0;
      model_traits.surv_resid[i] = 1.0 - (suscept * aff_cdf_term);
    }
    else // Not affected
    {
      model_traits.cond_probs[i] =   (suscept - suscept * unaff_cdf_term)
                                   / (1.0 - suscept * unaff_cdf_term);
      model_traits.surv_resid[i] = -suscept * unaff_cdf_term;
    }
  }
}

} // End namespace AO
} // End namespace SAGE


