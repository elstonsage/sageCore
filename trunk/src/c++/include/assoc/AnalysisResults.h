#ifndef ANALYSIS_RESULTS_H
#define ANALYSIS_RESULTS_H
///======================================================================
///
///  File:	AnalysisResults.h
///
///  Author:	Stephen Gross
///
///  Copyright 2002 R. C. Elston
///======================================================================

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include "maxfun/maxfun.h"
#include "numerics/cephes.h"
#include "numerics/fmatrix.h"
#include "maxfunapi/maxfunapi.h"
#include "assoc/OutputFormatter.h"
#include "containers/HashTable.h"
#include "sampling/sampling.h"
#include "assoc/Datatypes.h"
#include "assoc/Configuration.h"
#include "assoc/Residuals.h"

namespace SAGE  {
namespace ASSOC {

class AnalysisResults
{
  public:
  
    struct ModelResults
    {
      ModelResults();
    
      bool  sameVarianceComponents() const;
    
      string                                name;
      OUTPUT::Table                         sample_summary;
      vector<SAMPLING::Field::SummaryInfo>  covariate_infos;
      MAXFUN::Results                       maxfun_results_null;
      MAXFUN::Results                       maxfun_results_alt;
      MAXFUN::ParameterMgr                  score_test_parameter_mgr;
      double                                score_statistic;
      Residuals                             residuals_null;
      Residuals                             residuals_alt;
      Residuals                             ind_residuals_null;
      Residuals                             ind_residuals_alt;
      double                                affected_mean_residuals_null;
      double                                unaffected_mean_residuals_null;
      double                                affected_mean_residuals_alt;
      double                                unaffected_mean_residuals_alt;
      
      double                                residuals_scaling_factor;
    };
    
    // Summary information a model in which there is ONE test covariate.
    // The last covariate in the vector is the test/alternate covariate.
    struct ComparisonInfo
    {
      ComparisonInfo();
      double  waldPvalue() const;
      double  largerPvalue() const;
      double  pvalueRatio() const;
    
      string                     model_name;                // Name of the model
      bool                       converged;                 // Whether or not both NULL and ALT converged
      bool                       baseline_converged;        // Whether NULL converged
      MAXFUN::Parameter          intercept;                 // The final estimates for the intercept
      vector<MAXFUN::Parameter>  covs;                      // The final estimates for the covariate coefficients
      double                     lrt_pvalue;                // The p-value for the likelihood ratio test
      double                     score_pvalue;              // The p-value for the score test      
      bool                       same_variance_components;  // Do alt and null models have the same effects?
    };
    
    class LRTLess
    {
      public:
        bool  operator()(const ComparisonInfo& lhs, const ComparisonInfo& rhs) const;
    };
    
    class WaldLess
    {
      public:
        bool  operator()(const ComparisonInfo& lhs, const ComparisonInfo& rhs) const;    
    };
    
    class LargerLess
    {
      public:
        bool  operator()(const ComparisonInfo& lhs, const ComparisonInfo& rhs) const;    
    };
    
    class RatioLess
    {
      public:
        bool  operator()(const ComparisonInfo& lhs, const ComparisonInfo& rhs) const;    
    };    
    
    class LRTFilter
    {
      public:
        LRTFilter(double cutpoint);
        bool  operator()(const ComparisonInfo& c) const;
        
      private:
        double  my_cutpoint;
    };
    
    class WaldFilter
    {
      public:
        WaldFilter(double cutpoint);
        bool  operator()(const ComparisonInfo& c) const;
        
      private:
        double  my_cutpoint;
    };
    
    static bool  less(double lhs, double rhs);
    static bool  greater(double lhs, double rhs);
    
    typedef vector<ComparisonInfo> ComparisonInfoVector;
    typedef vector<OUTPUT::Section> SectionVector;
  
    AnalysisResults(const Configuration& config);
    
    void addModelResults(const ModelResults& r, const Sampledata& sampledata);
  
    // Accessing data
    const Configuration& getConfig() const;
    const map<string, string>&  getNullResidualOutputs() const;
    const map<string, string>&  getAltResidualOutputs() const;
    
    // Public utility functions:
    void  generateDetailedOutput(ofstream& out);
    void  generateSummaryOutput(ofstream& out);
    void  generate_cov_matrix(ofstream&) const;
    void  generateTsvOutput(ofstream& tsv_file);

  private:
    AnalysisResults(const AnalysisResults& other);
    AnalysisResults& operator=(const AnalysisResults& other);
  
    // Private utility functions
    string bool_out(bool) const;
    

    
    OUTPUT::Table summarizeConfig() const;
      void  summarizeDisplayOptions(OUTPUT::Table& options) const;
    
    void  addModelDetailedSection(const ModelResults& results);
    
    void  addModelResidualOutput(const ModelResults& results, const Sampledata& sampledata);
      void  write_residual(ostringstream& out, double residual, double ind_residual, const Sampledata& sampledata, size_t member_id) const;
      const string  determine_sex(const FPED::Member& member) const;
    
    /*
    void  addIndependentResiduals(const ModelResults& results, const Sampledata& sampledata);
      void  calculateIndependentResiduals(const MAXFUN::ParameterMgr& p_manager,
                const vector<vector<vector<double > > >& shared_effects,
                map<string, string>& indep_residual_results,
                const Sampledata& sampledata, const string& model_name);
    */
    void  addPValues(const ModelResults& results); 
    
    // Turns the ComparisonInfoVector into a summary table.
    OUTPUT::Table  generateSummaryTable(vector<string>& no_lrts, bool complete) const;
      void  removeInvalidResults(list<ComparisonInfo>& results, vector<string>& no_lrts) const;
      void  orderResults(list<ComparisonInfo>& results) const;
      void  filterResults(list<ComparisonInfo>& results) const;
    
    // Data members
    const Configuration&  my_config;
    ComparisonInfoVector  my_cinfos;
    SectionVector   my_detailed_sections;
    map<string, string>  my_null_residual_outputs;           // model name, residuals
    map<string, string>  my_alt_residual_outputs;            // model name, residuals
    map<string, string>  my_null_indep_residual_outputs;     // model name, residuals
    map<string, string>  my_alt_indep_residual_outputs;      // model name, residuals    
};

}
} 

#endif
