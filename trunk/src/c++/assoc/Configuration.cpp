#include "assoc/Configuration.h"

namespace SAGE  {
namespace ASSOC {

const double  VARIANCE_COMPONENT_LB = 0.0000001;

string
displayOrder2String(Configuration::DisplayOrderEnum order)
{
  string  order_string;

  switch(order)
  {
    case Configuration::AS_INPUT:
      order_string = "As input";
      break;
    case Configuration::BY_LRT:
      order_string = "By LRT p-value";
      break;
    case Configuration::BY_WALD:
      order_string = "By Wald p-value";
      break;
    case Configuration::BY_LARGER_PVALUE:
      order_string = "by larger of LRT and Wald p-values";
      break;
    case Configuration::BY_PVALUE_RATIO:
      order_string = "by ratio of LRT and Wald p-values";
      break;
    default:
      assert(false);
  }

  return  order_string;
}

Configuration::effect
string2effect(const string& effect_name)
{

  Configuration::effect converted_effect = Configuration::RANDOM;
  if(effect_name == "Random")
  {
    converted_effect = Configuration::RANDOM;
  }
  else if(effect_name == "Polygenic")
  {
    converted_effect = Configuration::POLYGENIC;
  }
  else if(effect_name == "Family")
  {
    converted_effect = Configuration::FAMILY;
  }
  else if(effect_name == "Sibling")
  {
    converted_effect = Configuration::SIBLING;
  }
  else if(effect_name == "Marital")
  {
    converted_effect = Configuration::MARITAL;
  }
  else
  {
    assert(false);
  }
  return(converted_effect);
}


//========================================================================
// IMPLEMENTATION: Configuration::Model
//========================================================================
//
void
Configuration::Model::dump() const
{
  cout << "\n" << name << endl;
  
  cout << boolalpha;
  cout << "Write residuals:  " << endl;
  cout << "null residuals   " << null_residuals << endl;
  cout << "test residuals    " << test_residuals << endl;
  cout << noboolalpha;

  cout << "Null variance components:  " << endl;
  map<string, bool>::const_iterator  null_iter     = null_effects.begin();
  map<string, bool>::const_iterator  null_end_iter = null_effects.end();
  for(; null_iter != null_end_iter; ++null_iter)
  {
    if(null_iter->second)
    {
      cout << null_iter->first << "  ";
    }
  }

  cout << endl;

  cout << "Alt variance components:   " << endl;
  map<string, bool>::const_iterator  alt_iter     = alt_effects.begin();
  map<string, bool>::const_iterator  alt_end_iter = alt_effects.end();
  for(; alt_iter != alt_end_iter; ++alt_iter)
  {
    if(alt_iter->second)
    {
      cout << alt_iter->first << "  ";
    }
  }

  cout << endl;

  cout << "Null indices:  ";
  for(size_t i = 0; i < null_covariate_indices.size(); ++i)
  {
    cout << null_covariate_indices[i] << "  ";
  }

  cout << "\nTest indices:  ";
  for(size_t i = 0; i < alt_covariate_indices.size(); ++i)
  {
    cout << alt_covariate_indices[i] << "  ";
  }

  cout << endl;
}


//========================================================================
// IMPLEMENTATION: Configuration
//========================================================================
//
Configuration::Configuration()
{
  reverse_sort = false;

  my_title                = "";
  my_primary_trait_name   = "";
  my_maxfun_debug         = false;
  my_ofilename            = "";
  my_allow_averaging      = false;
  my_omit_complete_summary = false;
  my_dependent_trait_type = QUANTITATIVE;

  my_effects.push_back(Parameter("Variance components", "Random",
                                 MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY));
  my_effects.push_back(Parameter("Variance components", "Polygenic",
                                 MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY));
  my_effects.push_back(Parameter("Variance components", "Family",
                                 MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY));
  my_effects.push_back(Parameter("Variance components", "Sibling",
                                 MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY));
  my_effects.push_back(Parameter("Variance components", "Marital",
                                 MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY));

  my_transf_config.set_type(MFSUBMODELS::Transformation::GEORGE_ELSTON);
  my_transf_config.set_lambda1_fixed(false);
  my_transf_config.set_lambda1_init_est(1.0);
  my_transf_config.set_lambda1_lower_bound(-numeric_limits<double>::infinity());
  my_transf_config.set_lambda1_upper_bound(numeric_limits<double>::infinity());
  my_transf_config.set_lambda1_pvalue_mean(1.0);
  
  my_transf_config.set_lambda2_fixed(true);
  my_transf_config.set_lambda2_init_est(0.0);
  my_transf_config.set_lambda2_pvalue_mean(1.0);

  models_initialized = false;
  my_models.insert(make_pair("Baseline", Model("Baseline")));
  my_model_list.push_back("Baseline");

  my_display_order           = BY_LRT;

  my_display_all             = false;
  my_filter_by_wald          = false;
  my_filter_by_lrt           = false;
  my_limit_number_displayed  = true;
  my_wald_filter_threshold   = .01;
  my_lrt_filter_threshold    = .01;
  my_display_number          = 10;

  my_last_fixed_effect = my_effects.back().param_name;
  transform_both_sides = false;
  score_test_only = false;
}

Configuration::Configuration(const Configuration& other)
{
  reverse_sort                = other.reverse_sort;
  my_title                    = other.my_title;
  my_primary_trait_name       = other.my_primary_trait_name;
  my_maxfun_debug             = other.my_maxfun_debug;
  my_ofilename                = other.my_ofilename;
  my_allow_averaging          = other.my_allow_averaging;
  my_dependent_trait_type     = other.my_dependent_trait_type;
  my_effects                  = other.my_effects;
  my_cov_infos                = other.my_cov_infos;
  my_omit_complete_summary    = other.my_omit_complete_summary;
  my_transf_config            = other.my_transf_config;
  models_initialized          = other.models_initialized;
  my_display_order            = other.my_display_order;
  my_display_all              = other.my_display_all;
  my_filter_by_wald           = other.my_filter_by_wald;
  my_filter_by_lrt            = other.my_filter_by_lrt;
  my_limit_number_displayed   = other.my_limit_number_displayed;
  my_wald_filter_threshold    = other.my_wald_filter_threshold;
  my_lrt_filter_threshold     = other.my_lrt_filter_threshold;
  my_display_number           = other.my_display_number;
  my_model_list               = other.my_model_list;
  my_models                   = other.my_models;
  my_init_null_covariates     = other.my_init_null_covariates;
  my_last_fixed_effect        = other.my_last_fixed_effect;
  transform_both_sides        = other.transform_both_sides;
  score_test_only             = other.score_test_only;
}

Configuration&
Configuration::operator=(const Configuration& other)
{
  if(this != &other)
  {
    reverse_sort                = other.reverse_sort;
    my_title                    = other.my_title;
    my_primary_trait_name       = other.my_primary_trait_name;
    my_maxfun_debug             = other.my_maxfun_debug;
    my_ofilename                = other.my_ofilename;
    my_allow_averaging          = other.my_allow_averaging;
    my_dependent_trait_type     = other.my_dependent_trait_type;
    my_effects                  = other.my_effects;
    my_cov_infos                = other.my_cov_infos;
    models_initialized          = other.models_initialized;
    my_omit_complete_summary    = other.my_omit_complete_summary;
    my_transf_config            = other.my_transf_config;
    my_display_order            = other.my_display_order;
    my_display_all              = other.my_display_all;
    my_filter_by_wald           = other.my_filter_by_wald;
    my_filter_by_lrt            = other.my_filter_by_lrt;
    my_limit_number_displayed   = other.my_limit_number_displayed;
    my_wald_filter_threshold    = other.my_wald_filter_threshold;
    my_lrt_filter_threshold     = other.my_lrt_filter_threshold;
    my_display_number           = other.my_display_number;
    my_model_list               = other.my_model_list;
    my_models                   = other.my_models;
    my_init_null_covariates     = other.my_init_null_covariates;
    my_last_fixed_effect        = other.my_last_fixed_effect;
    transform_both_sides        = other.transform_both_sides;
    score_test_only             = other.score_test_only;    
  }

  return *this;
}


const vector<string>&
Configuration::getModelList() const
{
  return  my_model_list;
}

void
Configuration::removeBaselineModel()
{
  my_model_list.erase(my_model_list.begin());
  my_models.erase(my_models.find("Baseline"));
}

bool
Configuration::hasCovariate(const string& name) const
{
  return  find_if(my_cov_infos.begin(), my_cov_infos.end(), isNamed(name)) != my_cov_infos.end();
}

bool
Configuration::hasNullCovariate(const string& name) const
{
  AllNullCovConstIter  cov_iter     = allNullCovariateBegin();
  AllNullCovConstIter  cov_end_iter = allNullCovariateEnd();
  for(; cov_iter != cov_end_iter; ++cov_iter)
  {
    if(isNamed(name)(*cov_iter))
    {
      return  true;
    }
  }

  return  false;
}

const Configuration::CovariateInfo&
Configuration::getCovariateInfo(const string& name) const
{
  CovariateInfoVector::const_iterator  end_iter = my_cov_infos.end();
  CovariateInfoVector::const_iterator  iter     = find_if(my_cov_infos.begin(), end_iter, isNamed(name));

  assert(iter != end_iter);

  return  *iter;
}

size_t
Configuration::getAltCovCount(const string& model_name) const
{
  assert(my_models.count(model_name));

  return  my_models.find(model_name)->second.alt_covariate_indices.size();
}

void
Configuration::addCovariate(const CovariateInfo& covariate, const ModelList& models)
{
  /*
  cout << "\n\n\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Adding covariate " << covariate.cfg.param_name << " ..." << endl;
  cout << "With models: ";
  for(size_t i = 0; i < models.size(); ++i)
  {
    cout << models[i] << "  ";
  }

  cout << endl;
  */

  my_cov_infos.push_back(covariate);
  size_t  covariate_index = my_cov_infos.size() - 1;

  if(models.empty())
  {
    my_init_null_covariates.push_back(covariate_index);
  }
  else
  {
    ModelList::const_iterator  name_iter     = models.begin();
    ModelList::const_iterator  name_end_iter = models.end();
    for(; name_iter != name_end_iter; ++name_iter)
    {
      string  model_name = *name_iter;
      map<string, Model>::iterator  model_iter = my_models.find(model_name);
      if(model_iter == my_models.end())
      {
        Model  new_model(model_name);
        new_model.alt_covariate_indices.push_back(covariate_index);
        my_models.insert(make_pair(model_name, new_model));
        my_model_list.push_back(model_name);
      }
      else
      {
        model_iter->second.alt_covariate_indices.push_back(covariate_index);
      }
    }
  }

  /*
  cout << "Finished adding covariate " << covariate.cfg.param_name << " ..." << endl;
  dumpCovariates();
  dumpModels();
  for(size_t i = 0; i < my_init_null_covariates.size(); ++i)
  {
    cout << my_init_null_covariates[i];
  }

  cout << endl;
  */
}

// - Populate models with null covariate indices and variance component inclusion
//   information.
//
void
Configuration::initializeModels(bool polygenic_eff, bool sibling_eff, bool marital_eff, bool family_eff)
{
  map<string, Model>::iterator  m_iter     = my_models.begin();
  map<string, Model>::iterator  m_end_iter = my_models.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    supplyVCStatus(m_iter, polygenic_eff, sibling_eff, marital_eff, family_eff);
    supplyNullIndices(m_iter);
  }

  models_initialized = true;
}

// - Populate model with the inclusion status of each variance component.
//
void
Configuration::supplyVCStatus(map<string, Model>::iterator& model_iter,
                              bool polygenic_eff, bool sibling_eff, bool marital_eff, bool family_eff)
{
  assert(! models_initialized);

  // - Fixed effects.
  model_iter->second.null_effects["Polygenic"] = polygenic_eff;
  model_iter->second.null_effects["Family"] = family_eff;
  model_iter->second.null_effects["Marital"] = marital_eff;
  model_iter->second.null_effects["Sibling"] = sibling_eff;

  model_iter->second.alt_effects["Polygenic"] = polygenic_eff;
  model_iter->second.alt_effects["Family"] = family_eff;
  model_iter->second.alt_effects["Marital"] = marital_eff;
  model_iter->second.alt_effects["Sibling"] = sibling_eff;

  // - User specified effects.
  size_t  total_effect_count = my_effects.size();
  size_t  ue = model_iter->second.null_effects.size() + 1;   // + 1 accounts for random effect which
                                                             // which is part of my_effects but not part
                                                             // model effects.
  for(; ue < total_effect_count; ++ue)
  {
    const string&  user_effect_name = my_effects[ue].param_name;
    model_iter->second.null_effects[user_effect_name] = true;
    model_iter->second.alt_effects[user_effect_name] = true;
  }
}

// - Populate model with indices of null covariates.
//
void
Configuration::supplyNullIndices(map<string, Model>::iterator& model_iter)
{
  model_iter->second.null_covariate_indices = my_init_null_covariates;
}

OUTPUT::Section
Configuration::dump() const
{
  OUTPUT::Section  s("Configuration");

  s << (OUTPUT::Table()
    << (OUTPUT::TableRow() << "title" << my_title)
    << (OUTPUT::TableRow() << "primary trait" << my_primary_trait_name)
    << (OUTPUT::TableRow() << "ofilename" << my_ofilename)
    << (OUTPUT::TableRow() << "allow averaging?" << my_allow_averaging)
    << (OUTPUT::TableRow() << "maxfun debug?" << my_maxfun_debug)
    << (OUTPUT::TableRow() << "dependent trait type" << (size_t)my_dependent_trait_type)
    << (OUTPUT::TableRow() << "display order" << displayOrder2String(my_display_order))
    << (OUTPUT::TableRow() << "display all?" << my_display_all)
    << (OUTPUT::TableRow() << "filter by wald?" << my_filter_by_wald)
    << (OUTPUT::TableRow() << "filter by LRT?" << my_filter_by_lrt)
    << (OUTPUT::TableRow() << "limit number displayed?" << my_limit_number_displayed)
    << (OUTPUT::TableRow() << "wald filter threshold" << my_wald_filter_threshold)
    << (OUTPUT::TableRow() << "lrt threshold" << my_lrt_filter_threshold)
    << (OUTPUT::TableRow() << "display number" << my_display_number));

  vector<string>  model_list = getModelList();

  OUTPUT::Table  t("Model list");

  for(std::vector<std::string>::const_iterator m = model_list.begin(); m != model_list.end(); ++m)
    t << (OUTPUT::TableRow() << *m);

  s << t;

  return s;
}

void
Configuration::dumpModels() const
{
  OUTPUT::Section  s("Models");

  map<string, Model>::const_iterator  m_iter     = my_models.begin();
  map<string, Model>::const_iterator  m_end_iter = my_models.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    m_iter->second.dump();
  }
}

void
Configuration::dumpCovariates() const
{
  for(size_t i = 0; i < my_cov_infos.size(); ++i)
  {
    cout << i << "  " << my_cov_infos[i].cfg.param_name << endl;
  }
}

}
}

