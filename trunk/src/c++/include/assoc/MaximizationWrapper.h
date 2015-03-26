#ifndef ASSOC_ANALYSIS_WRAPPER_H
#define ASSOC_ANALYSIS_WRAPPER_H
//==============================================================================
//
//  File:	MaximizationWrapper.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//==============================================================================


#include <math.h>
#include <cmath>
#include <fstream>
#include "numerics/fmatrix.h"
#include "numerics/sinfo.h"
#include "LSF/parse_ops.h"
#include "app/output_streams.h"
#include "fped/fped.h"
#include "maxfunapi/maxfunapi.h"
#include "maxfun/maxfun.h"
#include "maxfun/sub_model.h"
#include "assoc/AnalysisResults.h"
#include "assoc/Datatypes.h"
#include "assoc/Configuration.h"
#include "assoc/Calculator.h"

namespace SAGE  {
namespace ASSOC {

const double  BINARY_DIVISOR = 1.0;     // Determined emperically.

class vc_order
{
  public:
    vc_order();
    bool  operator()(const string& left, const string& right);
  
  private:
    map<string, int>  my_order;  
};

/// The analysis wrapper is bound to one specific analysis. It initiates that
/// analysis and returns the results to the calling object.
///
class MaximizationWrapper
{
  public:
    typedef MFSUBMODELS::Transformation::Configuration  TransConfig;
    static void maximize(const std::string& model_name, Configuration& config,
                         const FPED::Multipedigree& fped, const Sampledata& sampledata,
                         AnalysisResults::ModelResults& results, ostream& messages,
                         bool base_initialized, cerrorstream  errors = sage_cerr);

  private:
    MaximizationWrapper(const std::string& model_name, Configuration& config,
                        const FPED::Multipedigree& fped, const Sampledata& sampledata,
                        AnalysisResults::ModelResults& results, ostream& messages, bool base_initialized, 
                        cerrorstream  errors = sage_cerr);

    void  doMaximization(AnalysisResults::ModelResults& results, 
                          const TransConfig& trans_config);
      void  setMaximizationSequence(MAXFUN::SequenceCfg& cfg) const;
      void fixToNoTrans(TransConfig& trans_config) const;
      void  setUsersCovVals(MAXFUN::ParameterMgr& mgr) const;
      void  checkResults(const MAXFUN::Results& results) const;
      string  checkVarianceComponents(const MAXFUN::ParameterMgr& mgr);
        string  checkSpecialVarianceComponents(const MAXFUN::ParameterMgr& mgr);
      void  removeVarianceComponent(const string& component);
      void  checkRandomVariance(const MAXFUN::ParameterMgr& mgr);
      
    void  doScoreTest(AnalysisResults::ModelResults& results, 
                              const TransConfig& trans_config);
      
    /// Adds all the necessary parameters and so on to the ParameterMgr.
    void addParams(AnalysisResults::ModelResults& results, MAXFUN::Function& func);
      void  addFixedEffects(MAXFUN::Function& func,  
                            vector<size_t>& fixed_effect_indices,
                            const vector<bool>& fixed_effect_usage);
      void  addUserEffects(MAXFUN::Function& func);
      void  addTotalVariance(AnalysisResults::ModelResults& results, MAXFUN::Function& func,  
                             const vector<size_t>& fixed_effect_indices,
                             const vector<bool>& fixed_effect_usage, size_t& tv_idx) const;
      void  addHeritability(MAXFUN::Function& func,  
                            const vector<size_t>& fixed_effect_indices,
                            const vector<bool>& fixed_effect_usage, size_t tv_idx) const; 
      void  addVarianceProportions(MAXFUN::Function& func,  
                                   const vector<size_t>& fixed_effect_indices,
                                   const vector<bool>& fixed_effect_usage, size_t tv_idx) const;                             
      void  addNullCovCoefficients(MAXFUN::Function& func) const;
      void  addTestCovCoefficients(MAXFUN::Function& func) const;
      void  addFamilialCorrelations(MAXFUN::Function& func, 
                                    const vector<size_t>& fixed_effect_indices,
                                    const vector<bool>& fixed_effect_usage, size_t tv_idx) const;
      void  addEnvironmentalCorrelations(MAXFUN::Function& func, 
                                         const vector<size_t>& fixed_effect_indices,
                                         const vector<bool>& fixed_effect_usage, size_t tv_idx) const;
      void  addIntercept(MAXFUN::Function& func) const;      
    
    void  establishInitialEstimates(MAXFUN::ParameterMgr& mgr, double prevalence, 
                                    const AnalysisResults::ModelResults& results);
      void  calculateInitialEsts(MAXFUN::ParameterMgr& mgr, double prevalence);
        void  checkDataSufficiency(size_t num_of_valid_inds, const MAXFUN::ParameterMgr& mgr);
        void  populateMatrices(FortranMatrix<double>& main_phenotypes, FortranMatrix<double>& all_covariates, 
                               const MAXFUN::ParameterMgr& mgr, double prevalence);
          double  calculateBinaryPhenotype(double prevalence, size_t ind_idx, const MAXFUN::ParameterMgr& mgr) const;
            double  calculatePrevalence(size_t ind_idx, const MAXFUN::ParameterMgr& mgr) const;                             
        void  setVarianceComponents(double est_total_variance, MAXFUN::ParameterMgr& mgr) const;
          void  checkVarianceEstimates(MAXFUN::ParameterMgr& mgr) const;
      void  usePreviousResults(MAXFUN::ParameterMgr& mgr, const MAXFUN::ParameterMgr& previous_mgr);
        void  usePreviousLambdas(MAXFUN::ParameterMgr& mgr, const MAXFUN::ParameterMgr& previous_mgr) const;    
       
  
      void  rememberEstimates(vector<double>& estimates, const MAXFUN::ParameterMgr& mgr) const;
      void  setNewInitialValues(MAXFUN::ParameterMgr& mgr);
      double  calculateTotalBetaChange(const vector<double>& previous, const vector<double>& current) const;
      
    // Independent residuals
    Residuals  calculateIndResiduals(const MAXFUN::ParameterMgr& mgr, const Calculator& calc) const;
      void  populateVarCovMatrix(FortranMatrix<double>& var_cov,
                                 const map<string, FortranMatrix<double> >& shared_effects,
                                 const MAXFUN::ParameterMgr& mgr, size_t member_count      ) const;
      void populateExpValues(FortranMatrix<double>& exp_values, const vector<double>& diffs, 
                             const map<FPED::MemberConstPointer, size_t>& index_lookup, size_t member_count) const;
      void populatePhenotypes(FortranMatrix<double>& phenotypes, 
                               const map<FPED::MemberConstPointer, size_t>& index_lookup, size_t member_count) const;
      void reindexResiduals(vector<double>& ind_residuals_vector,
                            const FortranMatrix<double>& ind_residuals,
                            const vector<FPED::MemberConstPointer>& member_lookup) const;
                                                                                  
                                                                                                 
                                                                                                                                       
    
    // Data members
    string                      my_model_name;
    Configuration&              my_config;
    ostream&                    my_messages;
    bool                        null_results_available;
    const Sampledata&           my_sampledata;
    const FPED::Multipedigree&  my_fped;
    cerrorstream                my_errors;
    bool                        alt_covariates_included;
    bool                        score_test;
    double                      my_divisor;
    
    enum  max_phase { ONE_OF_TWO, TWO_OF_TWO, ONE_OF_ONE } my_phase;
};

}
} 

#endif
