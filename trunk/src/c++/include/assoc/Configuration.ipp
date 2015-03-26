//============================================================================
// File:      Configuration.ipp
//
// Author:    Stephen Gross
//
// History:   11/12/7  moved and cleaned up inline functions from 
//                     Configuration.h        -djb
//
// Notes:     
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================


//============================================================================
// IMPLEMENTATION:  Configuration::Model
//============================================================================
//


inline
Configuration::Model::Model(const string& model_name)
    : name(model_name), null_residuals(false), test_residuals(false)
{
  null_effects["Polygenic"] = true; 
  null_effects["Family"] = false;
  null_effects["Marital"] = true;
  null_effects["Sibling"] = true;
  alt_effects["Polygenic"] = true; 
  alt_effects["Family"] = false;
  alt_effects["Marital"] = true; 
  alt_effects["Sibling"] = true;
}


//============================================================================
// IMPLEMENTATION:  Configuration::NullCovConstIter
//============================================================================
//
inline
Configuration::NullCovConstIter::NullCovConstIter(const CovariateInfoVector& covariates, 
                                                            const Model& model, bool end)
    : my_covariates(&covariates), my_end_iterator(model.null_covariate_indices.end()),
      my_iterator(end ? my_end_iterator : model.null_covariate_indices.begin())
{}

inline
Configuration::NullCovConstIter::~NullCovConstIter()
{}

inline const Configuration::CovariateInfo&
Configuration::NullCovConstIter::operator *() const
{
  assert(my_iterator != my_end_iterator);
  assert(*my_iterator < my_covariates->size());
  
  return  (*my_covariates)[*my_iterator];
}

inline const Configuration::CovariateInfo*
Configuration::NullCovConstIter::operator ->() const
{
  return  &(**this);
}

inline Configuration::NullCovConstIter&
Configuration::NullCovConstIter::operator ++()
{
  ++my_iterator;
  
  return  *this;
}

inline bool
Configuration::NullCovConstIter::operator !=(const NullCovConstIter& other) const
{
  return  my_iterator != other.my_iterator       ||
          my_covariates != other.my_covariates;
}

//============================================================================
// IMPLEMENTATION:  Configuration::AllNullCovConstIter
//============================================================================
//
inline
Configuration::AllNullCovConstIter::AllNullCovConstIter(const CovariateInfoVector& covariates, 
                                                        const Configuration& config, bool end)
    : my_covariates(&covariates), my_end_iterator(config.my_init_null_covariates.end()),
      my_iterator(end ? my_end_iterator : config.my_init_null_covariates.begin())
{}

inline
Configuration::AllNullCovConstIter::~AllNullCovConstIter()
{}

inline const Configuration::CovariateInfo&
Configuration::AllNullCovConstIter::operator *() const
{
  assert(my_iterator != my_end_iterator);
  assert(*my_iterator < my_covariates->size());
  
  return  (*my_covariates)[*my_iterator];
}

inline const Configuration::CovariateInfo*
Configuration::AllNullCovConstIter::operator ->() const
{
  return  &(**this);
}

inline Configuration::AllNullCovConstIter&
Configuration::AllNullCovConstIter::operator ++()
{
  ++my_iterator;
  
  return  *this;
}

inline bool
Configuration::AllNullCovConstIter::operator !=(const AllNullCovConstIter& other) const
{
  return  my_iterator != other.my_iterator       ||
          my_covariates != other.my_covariates;
}

//============================================================================
// IMPLEMENTATION:  Configuration::TestCovConstIter
//============================================================================
//
inline
Configuration::TestCovConstIter::TestCovConstIter(const CovariateInfoVector& covariates, 
                                                            const Model& model, bool end)
    : my_covariates(&covariates), my_end_iterator(model.alt_covariate_indices.end()),
      my_iterator(end ? my_end_iterator : model.alt_covariate_indices.begin())
{}

inline
Configuration::TestCovConstIter::~TestCovConstIter()
{}

inline const Configuration::CovariateInfo&
Configuration::TestCovConstIter::operator *() const
{
  assert(my_iterator != my_end_iterator);
  assert(*my_iterator < my_covariates->size());
  
  return  (*my_covariates)[*my_iterator];
}

inline const Configuration::CovariateInfo*
Configuration::TestCovConstIter::operator ->() const
{
  return  &(**this);
}

inline Configuration::TestCovConstIter&
Configuration::TestCovConstIter::operator ++()
{
  ++my_iterator;
  
  return  *this;
}

inline bool
Configuration::TestCovConstIter::operator !=(const TestCovConstIter& other) const
{
  return  my_iterator != other.my_iterator       ||
          my_covariates != other.my_covariates;
}


//============================================================================
// IMPLEMENTATION:  Configuration::ModelCovConstIter
//============================================================================
//
inline
Configuration::ModelCovConstIter::ModelCovConstIter(const CovariateInfoVector& covariates, 
                                                            const Model& model, bool end)
    : my_covariates(&covariates), 
      my_end_iterator(model.alt_covariate_indices.end()),
      my_begin_alt_iterator(model.alt_covariate_indices.begin()),
      my_end_null_iterator(model.null_covariate_indices.end()),
      my_iterator(end ? my_end_iterator : model.null_covariate_indices.begin())
{
  if(my_iterator == my_end_null_iterator)
  {
    my_iterator = my_begin_alt_iterator;
  }
}

inline
Configuration::ModelCovConstIter::~ModelCovConstIter()
{}

inline const Configuration::CovariateInfo&
Configuration::ModelCovConstIter::operator *() const
{
  assert(my_iterator != my_end_iterator);
  assert(*my_iterator < my_covariates->size());
  
  return  (*my_covariates)[*my_iterator];
}

inline const Configuration::CovariateInfo*
Configuration::ModelCovConstIter::operator ->() const
{
  return  &(**this);
}

inline Configuration::ModelCovConstIter&
Configuration::ModelCovConstIter::operator ++()
{
  if(++my_iterator == my_end_null_iterator)
  {
    my_iterator = my_begin_alt_iterator;
  }
  
  return  *this;
}

inline bool
Configuration::ModelCovConstIter::operator !=(const ModelCovConstIter& other) const
{
  return  my_iterator != other.my_iterator       ||
          my_covariates != other.my_covariates;
}


//============================================================================
// IMPLEMENTATION:  Configuration::isNamed
//============================================================================
//
inline
Configuration::isNamed::isNamed(const string& name)
      : my_name(name)
{}

inline  bool  
Configuration::isNamed::operator()(const CovariateInfo& cov_info)  
{ 
  return  toUpper(my_name) == toUpper(cov_info.cfg.param_name); 
}


//============================================================================
// IMPLEMENTATION:  Configuration::parameterNamed
//============================================================================
//
inline
Configuration::parameterNamed::parameterNamed(const string& name)
      : my_name(name)
{}

inline  bool  
Configuration::parameterNamed::operator()(const Parameter& param)  
{ 
  return  toUpper(my_name) == toUpper(param.param_name); 
}


//============================================================================
// IMPLEMENTATION:  Configuration    
//============================================================================
//
inline void 
Configuration::setTitle(const std::string& title)               
{ 
  my_title = title; 
}

inline void 
Configuration::setPrimaryTraitName(const std::string& name)     
{ 
  my_primary_trait_name = name; 
}

inline void 
Configuration::setOfilename(const std::string& name)            
{ 
  my_ofilename = name; 
}
   
inline void 
Configuration::setAllowAveraging(bool a)                        
{ 
  my_allow_averaging = a; 
}

inline void 
Configuration::setDependentTraitType(DependentTraitTypeEnum d)  
{ 
  my_dependent_trait_type = d; 
}

inline void 
Configuration::setMaxfunDebug(bool d)                           
{ 
  my_maxfun_debug = d; 
}

inline void
Configuration::setTransformBothSides(bool both_sides)
{
  transform_both_sides = both_sides;
}

inline void
Configuration::setScoreTestOnly(bool score_only)
{
  score_test_only = score_only;
}

inline void 
Configuration::setOmitCompleteSummary(bool omit)                         
{ 
  my_omit_complete_summary = omit; 
}      

inline void 
Configuration::setTransformationOption(MFSUBMODELS::Transformation::TransformationType option)
{
  my_transf_config.set_type(option);
}

inline void  
Configuration::setDisplayOrder(DisplayOrderEnum order)
{
  my_display_order = order;
}

inline void  
Configuration::setDisplayAll(bool display_all)
{
  my_display_all = display_all;
}
        
inline void  
Configuration::setFilterByWald(bool filter_by_wald)
{
  my_filter_by_wald = filter_by_wald;
}
            
inline void  
Configuration::setFilterByLRT(bool filter_by_lrt)
{
  my_filter_by_lrt = filter_by_lrt;
}
                
inline void  
Configuration::setLimitNumberDisplayed(bool limit_number_displayed)
{
  my_limit_number_displayed = limit_number_displayed;
}
                    
inline void  
Configuration::setWaldFilterThreshold(double wald_filter_threshold)
{
  my_wald_filter_threshold = wald_filter_threshold;
}
                        
inline void  
Configuration::setLRTFilterThreshold(double lrt_filter_threshold)
{
  my_lrt_filter_threshold = lrt_filter_threshold;
}
                            
inline void  
Configuration::setDisplayNumber(size_t  display_number)
{
  my_display_number = display_number;
}

inline void
Configuration::removeEffectFromModel(const string& effect_name, const string& model_name, bool alt)
{
  assert(hasModel(model_name));
  assert(models_initialized);
  
  if(alt)
  {
    map<string, bool>&  alt_effects = my_models.find(model_name)->second.alt_effects;
    
    assert(alt_effects.count(effect_name) == 1);
    alt_effects[effect_name] = false;
  }
  else
  {
    map<string, bool>&  null_effects = my_models.find(model_name)->second.null_effects;
    
    assert(null_effects.count(effect_name) == 1);
    null_effects[effect_name] = false;  
  }
}

inline size_t
Configuration::getUserEffectCount(const string& model_name, bool alt) const
{
  assert(hasModel(model_name));
  assert(models_initialized);
  
  size_t  count = 0;
  
  const Configuration::Model&  model = my_models.find(model_name)->second;
  const map<string, bool>&  effects = alt ? model.alt_effects : model.null_effects;
  
  map<string, bool>::const_iterator  eff_iter     = effects.begin();
  map<string, bool>::const_iterator  eff_end_iter = effects.end();
  for(; eff_iter != eff_end_iter; ++eff_iter)
  {
    if(eff_iter->second)
    {
      ++count;
    }
  }
  return(count); // added per KCC to remove Solaris warning
}

inline void
Configuration::addUserEffect(const Parameter& eff)
{
  my_effects.push_back(eff);
}

inline Parameter&
Configuration::getFixedEffect(Configuration::effect eff)
{
  return  my_effects[eff];
}

inline Configuration::effect_iterator
Configuration::userEffectBegin()
{
  return  ++(find_if(my_effects.begin(), my_effects.end(), parameterNamed(my_last_fixed_effect)));
}

inline Configuration::effect_iterator
Configuration::userEffectEnd()
{
  return  my_effects.end();
}

inline Configuration::CovConstIter
Configuration::covariateBegin() const
{
  return  my_cov_infos.begin();
}

inline Configuration::CovConstIter
Configuration::covariateEnd() const
{
  return  my_cov_infos.end();
}

inline Configuration::NullCovConstIter
Configuration::nullCovariateBegin(const string& model_name) const
{
  assert(models_initialized);
  assert(my_models.count(model_name));

  return  NullCovConstIter(my_cov_infos, (my_models.find(model_name))->second, false);
}

inline Configuration::NullCovConstIter
Configuration::nullCovariateEnd(const string& model_name) const
{
  assert(models_initialized);
  assert(my_models.count(model_name));

  return  NullCovConstIter(my_cov_infos, (my_models.find(model_name))->second, true);
}

inline Configuration::AllNullCovConstIter
Configuration::allNullCovariateBegin() const
{
  return  AllNullCovConstIter(my_cov_infos, *this, false);
}

inline Configuration::AllNullCovConstIter
Configuration::allNullCovariateEnd() const
{
  return  AllNullCovConstIter(my_cov_infos, *this, true);
}

inline Configuration::TestCovConstIter
Configuration::testCovariateBegin(const string& model_name) const
{
  assert(my_models.count(model_name));

  return  TestCovConstIter(my_cov_infos, (my_models.find(model_name))->second, false);
}

inline Configuration::TestCovConstIter
Configuration::testCovariateEnd(const string& model_name) const
{
  assert(my_models.count(model_name));

  return  TestCovConstIter(my_cov_infos, (my_models.find(model_name))->second, true);
}


inline Configuration::ModelCovConstIter
Configuration::modelCovariateBegin(const string& model_name) const
{
  assert(models_initialized);
  assert(my_models.count(model_name));

  return  ModelCovConstIter(my_cov_infos, (my_models.find(model_name))->second, false);
}

inline Configuration::ModelCovConstIter
Configuration::modelCovariateEnd(const string& model_name) const
{
  assert(models_initialized);
  assert(my_models.count(model_name));

  return  ModelCovConstIter(my_cov_infos, (my_models.find(model_name))->second, true);
}


inline MFSUBMODELS::Transformation::Configuration&  
Configuration::getTransConfig()    
{ 
  return my_transf_config; 
}

inline const std::string& 
Configuration::getTitle() const                 
{ 
  return my_title; 
}

inline const std::string& 
Configuration::getPrimaryTraitName() const      
{ 
  return my_primary_trait_name; 
}

inline const std::string& 
Configuration::getOfilename() const             
{ 
  return my_ofilename; 
}

inline bool 
Configuration::getAllowAveraging() const                      
{ 
  return my_allow_averaging; 
}

inline bool 
Configuration::getMaxfunDebug() const 
{ 
  return my_maxfun_debug; 
}

inline bool
Configuration::getTransformBothSides() const
{
  return  transform_both_sides;
}

inline bool
Configuration::getScoreTestOnly() const
{
  return  score_test_only;
}

inline const DependentTraitTypeEnum& 
Configuration::getDependentTraitType() const                   
{ 
  return my_dependent_trait_type; 
}

inline const Parameter&
Configuration::getFixedEffect(Configuration::effect eff) const
{
  return  my_effects[eff];
}

inline Configuration::const_effect_iterator
Configuration::userEffectBegin() const
{
  return  ++(find_if(my_effects.begin(), my_effects.end(), parameterNamed(my_last_fixed_effect)));
}

inline Configuration::const_effect_iterator
Configuration::userEffectEnd() const
{
  return  my_effects.end();
}

inline bool
Configuration::modelHasFixedEffect(const string& model_name, bool alt, effect eff) const
{
  assert(models_initialized);
  map<string, Model>::const_iterator  model_iter = my_models.find(model_name);
  assert(model_iter != my_models.end());
  bool return_val = false;
  
  switch(eff)
  {
    case RANDOM:
      return_val = true;
      break;
    case POLYGENIC:
      return_val =  alt ? model_iter->second.alt_effects.find("Polygenic")->second :
                    model_iter->second.null_effects.find("Polygenic")->second;
      break;              
    case FAMILY:
      return_val =  alt ? model_iter->second.alt_effects.find("Family")->second :
                    model_iter->second.null_effects.find("Family")->second;      
      break;
    case SIBLING:
      return_val =  alt ? model_iter->second.alt_effects.find("Sibling")->second :
                    model_iter->second.null_effects.find("Sibling")->second;
      break;              
    case MARITAL:
      return_val =  alt ? model_iter->second.alt_effects.find("Marital")->second :
                    model_iter->second.null_effects.find("Marital")->second; 
      break;              
    default:
      assert(false);
  }
  
  return(return_val); // added per KCC to remove Solaris warning
}

inline bool
Configuration::modelHasEffect(const string& model_name, bool alt, const string& effect_name) const
{
  assert(models_initialized);
  map<string, Model>::const_iterator  model_iter = my_models.find(model_name);
  assert(model_iter != my_models.end());
  
  if(alt)
  {
    map<string, bool>::const_iterator  effect_iter = model_iter->second.alt_effects.find(effect_name);
    assert(effect_iter != model_iter->second.alt_effects.end());
    
    return  effect_iter->second;
  }
  else
  {
    map<string, bool>::const_iterator  effect_iter = model_iter->second.null_effects.find(effect_name);
    assert(effect_iter != model_iter->second.null_effects.end());
    
    return  effect_iter->second;  
  }
}

inline const MFSUBMODELS::Transformation::Configuration& 
Configuration::getTransConfig() const     
{ 
  return my_transf_config; 
}

inline bool 
Configuration::getOmitCompleteSummary() const 
{ 
  return my_omit_complete_summary; 
}

inline Configuration::DisplayOrderEnum
Configuration::getDisplayOrder() const  
{
  return  my_display_order;
}

inline bool              
Configuration::getDisplayAll() const
{
  return  my_display_all;
}

inline bool              
Configuration::getFilterByWald() const
{
  return  my_filter_by_wald;
}

inline bool              
Configuration::getFilterByLRT() const
{
  return  my_filter_by_lrt; 
}
 
inline bool              
Configuration::getLimitNumberDisplayed() const
{
  return  my_limit_number_displayed;
}
 
inline double            
Configuration::getWaldFilterThreshold() const
{
  return  my_wald_filter_threshold; 
}

inline double            
Configuration::getLRTFilterThreshold() const
{
  return  my_lrt_filter_threshold;  
}

inline size_t            
Configuration::getDisplayNumber() const
{
  return  my_display_number;
}

inline bool
Configuration::getModelsInitialized() const
{
  return  models_initialized;
}

inline bool
Configuration::hasModel(const string& model_name) const
{
  return  my_models.count(model_name);
}

inline Configuration::Model&
Configuration::getModel(const string& model_name)
{
  map<string, Model>::iterator  model_iter = my_models.find(model_name);
  assert(model_iter != my_models.end());
  
  return  model_iter->second;
}

inline const Configuration::Model&
Configuration::getModel(const string& model_name) const
{
  map<string, Model>::const_iterator  model_iter = my_models.find(model_name);
  assert(model_iter != my_models.end());
  
  return  model_iter->second;
}
