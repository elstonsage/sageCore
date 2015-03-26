//=======================================================================
//
//  File:  assoc.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include <stdlib.h>
#include <unistd.h>
#include <boost/functional.hpp>

#ifndef __WIN32__
  #include <sys/mman.h>
#endif

#include "assoc/assoc.h"

namespace SAGE  {
namespace ASSOC {


bool  
AssocIndividualValidator::isValid(size_t i, const SAMPLING::IndividualTraitData& trait_data) const
{
  SAMPLING::GroupInfoConstIterator group_info_itr      = trait_data.begin();
  SAMPLING::GroupInfoConstIterator group_info_end_itr  = trait_data.end();
  for(; group_info_itr != group_info_end_itr; ++group_info_itr)
  {
    string  group_name = group_info_itr->name;
    if(group_name != "Effects" && ! isValidGroup(trait_data, group_name))
       return  false;
  }
                                                
  return  true;
}

assoc::assoc(int _argc, char** _argv)
     : APP::SAGEapp(APP::APP_ASSOC, true, _argc, _argv)
{
  LSFInit();
}

int 
assoc::main()
{
  AppData data(name, debug());      
  print_title(data.info());        
  data.process_input(argc, argv);

  // - No one actually filtered here!
  //
  FPED::FilteredMultipedigree f(data.pedigrees());
  FPED::MPFilterer::add_multipedigree_filtered_by_subpedigrees(f, data.pedigrees(), FPED::always_keep());
  FPED::MPFilterer::add_multipedigree_filtered_by_unconnecteds(f, data.pedigrees(), FPED::always_keep());
  f.construct();

  if(f.pedigree_count() == 0)
  {
    data.errors() << priority(critical) << "Error: No valid data in sample! Exiting..."  << endl;

    exit(EXIT_FAILURE);
  }

  if(! data.getAnalyses().size())
  {
    data.errors() << priority(fatal) << "No valid analyses specified.  Program cannot continue." << endl;
    
    exit(EXIT_FAILURE);
  }
  
  perform_analyses(data, f);

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}


void
assoc::perform_analyses(AppData& data, const FPED::Multipedigree& fmp)
{
  // - For each analysis block in the parameter file
  //
  vector<Configuration>::iterator configuration_itr     = data.getAnalyses().begin();
  vector<Configuration>::iterator configuration_end_itr = data.getAnalyses().end();
  for(; configuration_itr != configuration_end_itr; ++configuration_itr)
  {
    bool  allow_averaging = configuration_itr->getAllowAveraging();
    MAXFUN::Results  baseline_results;
    bool  baseline_results_initialized = false;
  
    cout << "Performing analysis '" << configuration_itr->getTitle() << "'..." << endl;

    AnalysisResults  analysis_results(*configuration_itr);
    
    // - For each regression model
    //
    Configuration::ModelList  model_list = configuration_itr->getModelList();
    Configuration::ModelList::const_iterator  model_itr     = model_list.begin();
    Configuration::ModelList::const_iterator  model_end_itr = model_list.end();
    for(; model_itr != model_end_itr; ++model_itr)
    {
    
      cout << "  Setting up sample for model '" << *model_itr << "'..." << endl;
      
      Sampledata sampledata(fmp, data.errors());
        
      createSampleFields(sampledata, *configuration_itr, *model_itr);
      if(! checkSample(sampledata, fmp, *configuration_itr, *model_itr, data.errors()))
      {
        continue;
      }

      AnalysisResults::ModelResults  model_results;
      
      model_results.name            = *model_itr;
      model_results.sample_summary  =  createSampleSummary(sampledata);  // Altered 6-12-7. djb
      model_results.covariate_infos.clear();
      
      if(baseline_results_initialized)
      {
        model_results.maxfun_results_null = baseline_results;
        model_results.maxfun_results_null.setSequenceName(*model_itr + " without test covariates");
      }
  
      // - For each covariate in the regression model
      //
      SAMPLING::FieldConstIterator field_itr     = sampledata.getFieldBegin("Covariates");
      SAMPLING::FieldConstIterator field_end_itr = sampledata.getFieldEnd("Covariates");
      for(; field_itr != field_end_itr; ++field_itr)
      {
        model_results.covariate_infos.push_back(field_itr->getSummaryInfo());
      }
      
      //sampledata.dumpTraitValues();
      
      try
      {
        MaximizationWrapper::maximize(*model_itr, *configuration_itr, fmp, sampledata, 
                                      model_results, data.messages(), baseline_results_initialized, data.errors());
      }
      catch(const InsufficientData& e)
      {
        data.errors() << priority(error) << "Insufficient data to estimate parameters for model '" 
                      << e.what() << "'.  Try including 'allow_averaging=mean' in the analysis block "
                      << "of your parameter file or estimating fewer parameters.  Skipping model ..." << endl; 
                                  
        continue;        
      }
      catch(const BadLikelihood& e)
      {
        data.errors() << priority(error) << "Likelihood for model '" << e.what() << "' not "
                      << "finite.  Check for outliers and/or reduce the number of parameters "
                      << "to be estimated.  Skipping model ..." << endl;
                                  
        continue;
      }
      catch(const NonConvergenceWithTransformation& e)
      {
        data.errors() << priority(error) << "Maximization did not converge for model '" 
                      << e.what() << "'.  Skipping model ..." << endl; 
                                  
        continue;        
      }
     
      if(allow_averaging && (! baseline_results_initialized) && model_results.maxfun_results_null.getConverged())
      {
        baseline_results = model_results.maxfun_results_null;
        baseline_results_initialized = true;
      }
      
      analysis_results.addModelResults(model_results, sampledata);      
    }

    // Write output.
    ofstream  summary_ofile;
    ofstream  detailed_ofile;

    summary_ofile.open(string(configuration_itr->getOfilename() + ".sum").c_str());
    detailed_ofile.open(string(configuration_itr->getOfilename() + ".det").c_str());

    summary_ofile << getReleaseString();
    detailed_ofile << getReleaseString();

    analysis_results.generateSummaryOutput(summary_ofile);
    analysis_results.generateDetailedOutput(detailed_ofile);
    
    generateResidualOutput(*configuration_itr, analysis_results.getNullResidualOutputs(), NULL_RESIDUALS);
    generateResidualOutput(*configuration_itr, analysis_results.getAltResidualOutputs(), TEST_RESIDUALS);

    summary_ofile.close();
    detailed_ofile.close();
    
    #define CSV_OUTPUT
    #ifdef CSV_OUTPUT
      ofstream  tsv_file;
      
      tsv_file.open(string(analysis_results.getConfig().getOfilename() + ".tsv").c_str());
      analysis_results.generateTsvOutput(tsv_file);
      tsv_file.close();
    #endif
  } 
}


void
assoc::generateResidualOutput(const Configuration& config, const map<string, string>& residuals, residual_type r_type)
{
  string  base_filename = config.getOfilename();

  map<string, string>::const_iterator  r_iter     = residuals.begin();
  map<string, string>::const_iterator  r_end_iter = residuals.end();
  for(; r_iter != r_end_iter; ++r_iter)
  {
    const Configuration::Model&  model = config.getModel(r_iter->first);
    switch(r_type)
    {
      case NULL_RESIDUALS:
      {
        if(model.null_residuals)
        {
          ofstream  residuals_ofile;
          string  filename = base_filename + "_" + model.name + "_null.res";
          
          residuals_ofile.open(filename.c_str());
          residuals_ofile << r_iter->second;          
          residuals_ofile.close();
        }
        
        break;
      }
            
      case TEST_RESIDUALS:
      {
        if(model.test_residuals)
        {
          ofstream  residuals_ofile;
          string  filename = base_filename + "_" + model.name + "_test.res";
          
          residuals_ofile.open(filename.c_str());
          residuals_ofile << r_iter->second;
          residuals_ofile.close();
        }
          
        break;
      }
            
      default:
        assert(false);
    }
  }
}


void  
assoc::createSampleFields(SAMPLING::MemberDataSample& sample, const Configuration& config, const string& model_name) const
{
  sample.addGroup("Core traits");
  sample.addGroup("Covariates");
  sample.addGroup("Effects");

  // Import main phenotype
  //
  unsigned long  pheno_flags = SAMPLING::Field::NO_FLAGS;
  //unsigned long  pheno_flags = SAMPLING::Field::STDEV_ADJUST;
  sample.createField("Core traits", config.getPrimaryTraitName(), pheno_flags);
  sample.importField(config.getPrimaryTraitName(), "Core traits", "Main phenotype", pheno_flags);

  // Set up import flags for covariates
  // - To be valid individuals must have a value for the primary trait.  They must also have a value for every
  //   covariate OR allow_averaging must be specified by the user.  6-14-7 djb
  //
  //
  unsigned long import_flags = (config.getAllowAveraging() ? SAMPLING::Field::ALLOW_AVERAGING : SAMPLING::Field::NO_FLAGS) | 
                                SAMPLING::Field::STDEV_ADJUST                                                              | 
                                SAMPLING::Field::MEAN_ADJUST;

  // CREATE covariates.
  // 
  Configuration::ModelCovConstIter  cov_iter     = config.modelCovariateBegin(model_name);
  Configuration::ModelCovConstIter  cov_end_iter = config.modelCovariateEnd(model_name);
  for(; cov_iter != cov_end_iter; ++cov_iter)
  {
    sample.createField("Covariates", cov_iter->cfg.param_name, import_flags);
    sample.getField("Covariates", cov_iter->cfg.param_name).setUserDefined(false);
  }

  // POPULATE covariates
  //
  cov_iter = config.modelCovariateBegin(model_name);
  for(; cov_iter != cov_end_iter; ++cov_iter)
  {
      sample.importField(cov_iter->cfg.param_name, "Covariates", cov_iter->cfg.param_name, import_flags);
  }
  
  // CREATE and POPULATE user-specified effects. Note: class effects are user defined, 
  // but not the traits used to specify them.
  // 
  Configuration::const_effect_iterator  eff_iter     = config.userEffectBegin();
  Configuration::const_effect_iterator  eff_end_iter = config.userEffectEnd();
  for(; eff_iter != eff_end_iter; ++eff_iter)
  {
    sample.createField("Effects", eff_iter->param_name, SAMPLING::Field::NO_FLAGS);
    sample.getField("Effects", eff_iter->param_name).setUserDefined(false);
    sample.importField(eff_iter->param_name, "Effects", eff_iter->param_name, SAMPLING::Field::NO_FLAGS);    
  }  

  sample.finalizeData(AssocIndividualValidator());
  sample.finalizeUserCreatedData();
}

bool  
assoc::checkSample(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                   const Configuration& config, const string& model_name, cerrorstream& errors) const
{
  return  checkPrimaryTraitVariance(sample, model_name, config.getPrimaryTraitName(), errors) &&
          checkCovariateVariances(sample, model_name, errors)                                 && 
          checkForFamilies(sample, fped, config, model_name, errors)                          &&   
          checkForSibships(sample, fped, config, model_name, errors)                          && 
          checkForMatePairs(sample, fped, config, model_name, errors);                            
}
            

// - Make sure all covariates have non-zero variances.
//
bool
assoc::checkCovariateVariances(const SAMPLING::MemberDataSample& sample, 
                               const string& model_name, cerrorstream& errors) const
{
  SAMPLING::FieldConstIterator field_itr     = sample.getFieldBegin("Covariates");
  SAMPLING::FieldConstIterator field_end_itr = sample.getFieldEnd("Covariates");

  for(; field_itr != field_end_itr; ++field_itr)
  {
    //cout << "\nCovariate " << field_itr->getFieldName() << " variance  " << field_itr->getVariance() << endl;
  
    if(SAGE::isnan(field_itr->getStdev()) || field_itr->getStdev() < EFFECTIVELY_ZERO)
    {
      errors << priority(error) << "Covariate '" << field_itr->getFieldName() 
             << "' has no variance.  Skipping model '" << model_name << "' ..." << endl;
             
      return  false;
    }
  }
  
  return  true;
}


// - Make sure the primary trait has a non-zero variances.
//
bool
assoc::checkPrimaryTraitVariance(const SAMPLING::MemberDataSample& sample, const string& model_name, 
                                  const string& trait_name, cerrorstream& errors) const
{
  const SAMPLING::Field&  primary_trait = sample.getField("Core traits", "Main phenotype");
  
  //cout << "\nPramary trait variance  " << primary_trait.getVariance() << endl;

  if(SAGE::isnan(primary_trait.getStdev()) || primary_trait.getStdev() < EFFECTIVELY_ZERO)
  {
    errors << priority(error) << "Primary trait '" << trait_name 
           << "' has no variance for usable data for model '" << model_name 
           << "'.  Skipping model ..." << endl;
           
    return  false;
  }
  else
  {
    return  true;
  }
}


// If family effect enabled, verify that there is at least one family with two or more valid individuals.
//
bool  
assoc::checkForFamilies(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                        const Configuration& config, const string& model_name, cerrorstream& errors) const
{
  bool  valid = false;
  
  if(config.modelHasFixedEffect(model_name, true, Configuration::FAMILY) || config.modelHasFixedEffect(model_name, false, Configuration::FAMILY))
  {
    size_t  valid_fam_count = 0;

    // Loop across pedigrees
    for(FPED::PedigreeConstIterator ped_itr = fped.pedigree_begin(); ped_itr != fped.pedigree_end(); ++ped_itr)
    {
      if(valid)
      {
        break;
      }
          
      // Loop across families:
      for(FPED::FamilyConstIterator fam_itr = ped_itr->family_begin(); fam_itr != ped_itr->family_end(); ++fam_itr)
      {
      
        // Make a list of all the valid members
        size_t  valid_ind_count = sample.isValid(fam_itr->parent1()->mpindex()) +
                                  sample.isValid(fam_itr->parent2()->mpindex());

        for(FPED::OffspringConstIterator offspr_itr = fam_itr->offspring_begin(); offspr_itr != fam_itr->offspring_end(); ++offspr_itr)
        {
          valid_ind_count += sample.isValid(offspr_itr->mpindex());
        }
        
        valid_fam_count += valid_ind_count > 1;
        if(valid_fam_count > 1)
        {
          valid = true;
          break;
        }
      }
    }

    if(! valid)
    {
      errors << priority(error)
             << "This analysis includes family effect, but there are fewer than two informative families (that is, families"
             << " with at least two informative members. Skipping model '" << model_name << "' ..." << endl;
    }
  }
  else
  {
    valid = true;
  }
  
  return  valid;
}
                                       

// If sibling effect enabled, verify that there is at least one sibship with two or more valid individuals
//
bool  
assoc::checkForSibships(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,
                        const Configuration& config, const string& model_name, cerrorstream& errors) const
{
  bool  valid = false;

  if(config.modelHasFixedEffect(model_name, true, Configuration::SIBLING) || config.modelHasFixedEffect(model_name, false, Configuration::SIBLING))
  {
    // Loop across pedigrees
    for(FPED::PedigreeConstIterator ped_itr = fped.pedigree_begin(); ped_itr != fped.pedigree_end(); ++ped_itr)
    {
      if(valid)
      {
        break;
      }
    
      // Loop across families
      for(FPED::FamilyConstIterator fam_itr = ped_itr->family_begin(); fam_itr != ped_itr->family_end(); ++fam_itr)
      {
        // Count valid members.
        size_t  valid_count = 0;
        
        for(FPED::OffspringConstIterator offspr_itr = fam_itr->offspring_begin(); offspr_itr != fam_itr->offspring_end(); ++offspr_itr)
        {
          valid_count += sample.isValid(offspr_itr->mpindex());
        }
        
        if(valid_count > 1)
        {
          valid = true;
          break;
        }
      }
    }
    
    if(! valid)
    {
      errors << priority(error)
             << "This analysis includes sibling effect, but there are no sibships with at "
             << "least two informative members.  Skipping model '" << model_name << "' ..." << endl;
    }
  }
  else
  {
    valid = true; 
  }
  
  return  valid;
}


// If marital effect enabled, verify that there is at least one pair of mates in which both mates are valid.
//
bool  
assoc::checkForMatePairs(const SAMPLING::MemberDataSample& sample, const FPED::Multipedigree& fped,  
                         const Configuration& config, const string& model_name, cerrorstream& errors) const
{
  bool  valid = false;

  if(config.modelHasFixedEffect(model_name, true, Configuration::MARITAL) || config.modelHasFixedEffect(model_name, false, Configuration::MARITAL))
  {
    // Loop across pedigrees.
    for(FPED::PedigreeConstIterator ped_itr = fped.pedigree_begin(); ped_itr != fped.pedigree_end(); ++ped_itr)
    {
      if(valid)
      {
        break;
      }
    
      for(FPED::FamilyConstIterator fam_itr = ped_itr->family_begin(); fam_itr != ped_itr->family_end(); ++fam_itr)
      {
        if(sample.isValid(fam_itr->parent1()->mpindex()) && sample.isValid(fam_itr->parent2()->mpindex()))
        {
          valid = true;
          break;
        }
      }
    }
    
    if(! valid)
    {
      errors << priority(error)
             << "This analysis includes marital effect, but there are no families with two "
             << "informative parents.  Skipping model '" << model_name << "' ..." << endl;
    }
  }
  else
  {
    valid = true;
  }
  
  return  valid;
}
                                       

// - Added 6-12-7. djb
//
SAGE::OUTPUT::Table
assoc::createSampleSummary(const Sampledata& sd) const
{
  size_t  valid_singleton_count = sd.getValidSingletonCount();
  size_t  total_singleton_count = valid_singleton_count + sd.getInvalidSingletonCount();
  size_t  valid_subpedigree_member_count = sd.getValidIndividualCount() - valid_singleton_count;

  SAGE::OUTPUT::Table  sample_description("Sample description");
  sample_description << (SAGE::OUTPUT::TableRow() << "Number of individuals in dataset" << sd.getTotalIndividualCount())
                     << (SAGE::OUTPUT::TableRow() << "Number of constituent pedigrees in dataset" << getSubpedigreeCount(sd))
                     << (SAGE::OUTPUT::TableRow() << "Number of singletons in dataset" << total_singleton_count)
                     <<  SAGE::OUTPUT::Table::INSERT_BLANK_ROW()
                     << (SAGE::OUTPUT::TableRow() << "Number of constituent pedigree members with complete information" 
                                                  << valid_subpedigree_member_count)
                     << (SAGE::OUTPUT::TableRow() << "Number of singletons with complete information" << valid_singleton_count)
                     <<  SAGE::OUTPUT::Table::INSERT_BLANK_ROW();
                     
  return  sample_description;                    
}

// - Added 6-12-7. djb
//
size_t
assoc::getSubpedigreeCount(const Sampledata& sd) const
{
  size_t  count = 0;

  SAGE::FPED::PedigreeConstIterator  p_iter     = sd.getMultipedigree().pedigree_begin();
  SAGE::FPED::PedigreeConstIterator  p_end_iter = sd.getMultipedigree().pedigree_end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    count += p_iter->subpedigree_count();
  }
  
  return  count;
}

}
}

int main(int argc, char* argv[])
{
  SAGE::ASSOC::assoc assoc_instance(argc, argv);

  //lint -e(534) Ignoring return value
  assoc_instance.main();

  return 0;
}

