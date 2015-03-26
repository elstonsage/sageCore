#ifndef ASSOC_CONFIGURATION_H
#define ASSOC_CONFIGURATION_H
//=======================================================================
// File:	   Configuration.h
//
// Author:	 Stephen Gross
//
// History:  11/12/7  moved inlines to Configuration.ipp.   -djb
//           11/13/7  added ordering and filtering of displayed 
//                    results.  -djb
//           1/3/8    refactored to make iterating through cov-
//                    ariates more efficient.  Also made "model centric"
//                    rather than "covariate centric" for future flexibility.  -djb
//           3/27/8   Modified to accomodate user defined random effects.  -djb
//
// Copyright 2002 R. C. Elston
// All Rights Reserved
//=======================================================================


#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/cast.hpp>
#include "output/Output.h"
#include "maxfunapi/maxfunapi.h"
#include "mfsubmodels/mfsubmodels.h"
#include "assoc/Datatypes.h"

namespace SAGE  {
namespace ASSOC {

extern const double  VARIANCE_COMPONENT_LB;

typedef MAXFUN::ParameterInput  Parameter;

// - Describes the configuration options for a single ASSOC analysis.
//
class Configuration
{
  public:
    enum  effect { RANDOM, POLYGENIC, FAMILY, SIBLING, MARITAL };
  
    typedef vector<string> ModelList;

    struct CovariateInfo
    {
      Parameter  cfg;             // All the basic parameter options
    };
    
    typedef  vector<CovariateInfo> CovariateInfoVector;
    
    struct Model
    {
      Model(const string& model_name);
      void dump() const;
    
      string  name;
      bool  null_residuals;
      bool  test_residuals;
      vector<size_t>  null_covariate_indices;    
      vector<size_t>  alt_covariate_indices;
      
      // Random effects
      map<string, bool>  null_effects;
      map<string, bool>  alt_effects;
    };
    
    typedef CovariateInfoVector::const_iterator  CovConstIter;
    
    class NullCovConstIter  // Model specific
    {
      public:
        NullCovConstIter(const CovariateInfoVector& covariates, const Model& model, bool end);
        ~NullCovConstIter();
      
        const CovariateInfo&  operator  *() const;
        const CovariateInfo*  operator ->() const;
        NullCovConstIter&  operator ++();
        bool  operator !=(const NullCovConstIter& other) const;
        
      private:
        const CovariateInfoVector*  my_covariates;
        vector<size_t>::const_iterator  my_end_iterator;
        vector<size_t>::const_iterator  my_iterator;
    };
    
    class AllNullCovConstIter  // All null covariates for analysis block
    {
      public:
        AllNullCovConstIter(const CovariateInfoVector& covariates, 
                            const Configuration& config, bool end);
        ~AllNullCovConstIter();
      
        const CovariateInfo&  operator  *() const;
        const CovariateInfo*  operator ->() const;
        AllNullCovConstIter&  operator ++();
        bool  operator !=(const AllNullCovConstIter& other) const;
        
      private:
        const CovariateInfoVector*  my_covariates;
        vector<size_t>::const_iterator  my_end_iterator;
        vector<size_t>::const_iterator  my_iterator;
    };
    
    class TestCovConstIter 
    {
      public:
        TestCovConstIter(const CovariateInfoVector& covariates, const Model& model, bool end);
        ~TestCovConstIter();
      
        const CovariateInfo&  operator  *() const;
        const CovariateInfo*  operator ->() const;
        TestCovConstIter&  operator ++();
        bool  operator !=(const TestCovConstIter& other) const;
        
       private:
        const CovariateInfoVector*  my_covariates;
        vector<size_t>::const_iterator  my_end_iterator;
        vector<size_t>::const_iterator  my_iterator;
    };
    
    class ModelCovConstIter 
    {
      public:
        ModelCovConstIter(const CovariateInfoVector& covariates, const Model& model, bool end);
        ~ModelCovConstIter();
      
        const CovariateInfo&  operator  *() const;
        const CovariateInfo*  operator ->() const;
        ModelCovConstIter&  operator ++();
        bool  operator !=(const ModelCovConstIter& other) const;
        
       private:
        const CovariateInfoVector*  my_covariates;
        vector<size_t>::const_iterator  my_end_iterator;
        vector<size_t>::const_iterator  my_begin_alt_iterator;
        vector<size_t>::const_iterator  my_end_null_iterator;
        vector<size_t>::const_iterator  my_iterator;
    };    
    
    
    // - Added 10/3/7.   -djb
    class isNamed
    {
      public:
        isNamed(const string& name);
        
        bool  operator ()(const CovariateInfo& cov_info);
      
      private:
        string  my_name;
    };
    
    class parameterNamed
    {
      public:
        parameterNamed(const string& name);
        
        bool  operator ()(const Parameter& param);
      
      private:
        string  my_name;
    };    
    
    enum  DisplayOrderEnum  { AS_INPUT, BY_LRT, BY_WALD, BY_LARGER_PVALUE, BY_PVALUE_RATIO };
    typedef  vector<Parameter>::iterator  effect_iterator;
    typedef  vector<Parameter>::const_iterator  const_effect_iterator;        
    
    Configuration();
    Configuration(const Configuration& other);
    Configuration& operator=(const Configuration& other);
    
    // COVARIATES
    void  addCovariate(const CovariateInfo& covariate, const ModelList& models);
    
    CovConstIter  covariateBegin() const;         // All covariates in an analysis block
    CovConstIter  covariateEnd() const;           
    NullCovConstIter  nullCovariateBegin(const string& model_name) const;     // Null covariates in model
    NullCovConstIter  nullCovariateEnd(const string& model_name) const;       
    AllNullCovConstIter  allNullCovariateBegin() const;                       // All null covariates in analysis block
    AllNullCovConstIter  allNullCovariateEnd() const;           
    TestCovConstIter  testCovariateBegin(const string& model_name) const;     // Test covariates in model
    TestCovConstIter  testCovariateEnd(const string& model_name) const;
    ModelCovConstIter  modelCovariateBegin(const string& model_name) const;   // All covariates in model
    ModelCovConstIter  modelCovariateEnd(const string& model_name) const;
    bool  hasCovariate(const string& name) const;
    bool  hasNullCovariate(const string& name) const; 
    const CovariateInfo&  getCovariateInfo(const string& name) const;
    size_t  getAltCovCount(const string& model_name) const;                       
    
    // RANDOM EFFECTS
    void  addUserEffect(const Parameter& eff);
    Parameter&  getFixedEffect(effect eff);        
    void  removeEffectFromModel(const string& effect_name, const string& model_name, bool alt);
    
    size_t  getUserEffectCount(const string& model_name, bool alt) const;
    effect_iterator  userEffectBegin();
    effect_iterator  userEffectEnd();    
    const_effect_iterator  userEffectBegin() const;    
    const_effect_iterator  userEffectEnd() const;    
    const Parameter&  getFixedEffect(effect eff) const;    
    bool  modelHasFixedEffect(const string& model_name, bool alt, effect eff) const;    
    bool  modelHasEffect(const string& model_name, bool alt, const string& effect_name) const;        
    
    // MODELS
    bool  hasModel(const string& model_name) const;    
    const vector<string>&  getModelList() const;
    Model&  getModel(const string& model_name);
    const Model&  getModel(const string& model_name) const;
    void  removeBaselineModel();    
    
    // OUTPUT OPTIONS
    void  setOfilename(const string& name);    
    void  setTitle(const string& title);    
    void  setOmitCompleteSummary(bool omit);
    void  setDisplayOrder(DisplayOrderEnum order);
    void  setDisplayAll(bool display_all);
    void  setFilterByWald(bool filter_by_wald);
    void  setFilterByLRT(bool filter_by_lrt);
    void  setLimitNumberDisplayed(bool limit_number_displayed);
    void  setWaldFilterThreshold(double wald_filter_threshold);
    void  setLRTFilterThreshold(double lrt_filter_threshold);
    void  setDisplayNumber(size_t display_number);
    bool  getOmitCompleteSummary() const;
    DisplayOrderEnum  getDisplayOrder() const;
    bool  getDisplayAll() const;
    bool  getFilterByWald() const;
    bool  getFilterByLRT() const;
    bool  getLimitNumberDisplayed() const;
    double  getWaldFilterThreshold() const;
    double  getLRTFilterThreshold() const;
    size_t  getDisplayNumber() const;            
    
    const string&  getOfilename() const;
    const string&  getTitle() const;    
    
    // OTHER
    void  initializeModels(bool polygenic_eff, bool sibling_eff, bool marital_eff, bool family_eff);        
    void  setPrimaryTraitName(const string& name);
    void  setDependentTraitType(DependentTraitTypeEnum d);    
    void  setAllowAveraging(bool a);
    void  setMaxfunDebug(bool d);
    void  setTransformationOption(MFSUBMODELS::Transformation::TransformationType option);
    void  setTransformBothSides(bool both_sides);
    void  setScoreTestOnly(bool score_only);
    
    bool  getModelsInitialized() const;
    MFSUBMODELS::Transformation::Configuration&  getTransConfig();
    const MFSUBMODELS::Transformation::Configuration&   getTransConfig() const;    
    const DependentTraitTypeEnum&  getDependentTraitType() const;    
    const string&  getPrimaryTraitName() const;
    bool  getAllowAveraging() const;
    bool  getMaxfunDebug() const;
    bool  getTransformBothSides() const;
    bool  getScoreTestOnly() const;

    OUTPUT::Section  dump() const;
    void  dumpModels() const;
    void  dumpCovariates() const;
  
    // Debugging only! True means the reorganize_pds will simply reverse the order.
    bool reverse_sort;

  private:
    void  supplyVCStatus(map<string, Model>::iterator& model_iter,
                         bool polygenic_eff, bool sibling_eff, bool marital_eff, bool family_eff);
    void  supplyNullIndices(map<string, Model>::iterator& model_iter);  
  
    // Data members
    string  my_title;
    string  my_primary_trait_name;
    bool  my_maxfun_debug;
    string  my_ofilename;
    bool  my_allow_averaging;
    DependentTraitTypeEnum  my_dependent_trait_type;
    
    vector<Parameter>  my_effects;
    
    CovariateInfoVector  my_cov_infos;
    bool  my_omit_complete_summary;
    
    MFSUBMODELS::Transformation::Configuration  my_transf_config;
    
    // Trac #1764.  Summary display
    DisplayOrderEnum  my_display_order;
    
    bool    my_display_all;
    bool    my_filter_by_wald;
    bool    my_filter_by_lrt;
    bool    my_limit_number_displayed;
    double  my_wald_filter_threshold;
    double  my_lrt_filter_threshold;
    size_t  my_display_number;
    
    vector<string>  my_model_list;
    map<string, Model>  my_models;
    vector<size_t>  my_init_null_covariates;
    bool  models_initialized;
    string  my_last_fixed_effect;
    
    bool  transform_both_sides;
    bool  score_test_only;
};

string  displayOrder2String(Configuration::DisplayOrderEnum order);
Configuration::effect  string2effect(const string& effect_name);

#include "assoc/Configuration.ipp"

}
} 

#endif

