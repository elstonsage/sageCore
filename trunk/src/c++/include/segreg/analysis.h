#ifndef SEGREG_ANALYSIS_H
#define SEGREG_ANALYSIS_H

#include "segreg/segreg_calculator.h"
#include "segreg/PedigreeDataSet.h"
#include "segreg/model.h"
#include "app/output_streams.h"
#include "error/errorstream.h"
#include "boost/smart_ptr.hpp"

namespace SAGE
{

namespace SEGREG
{

class primary_analysis_results
{
    friend class primary_analysis;

  public:

    primary_analysis_results();
    primary_analysis_results(const primary_analysis_results&);

    primary_analysis_results& operator=(const primary_analysis_results&);

    ~primary_analysis_results();

    bool                     is_valid() const;

    const model&             get_final_model()         const;
    const MAXFUN::Results&   get_maxfun_results()      const;

    const FPED::Multipedigree&  get_multipedigree()       const;

    uint                     get_subpedigree_count()   const;
    uint                     get_unconnected_count()   const;

    /// If type probabilities are called for, we get them out with this
    /// function
    const pg::post_geno_map& get_type_probs() const;

    /// If pen_function probabilities are called for, we get them out with
    /// this function
    const pf::pen_func_map& get_pfunc_probs() const;

    const MlmResidCorrelationCalculator&     get_mlm_resid_corr() const;

    double get_one_mean_results(); // due to JA
    double get_one_susc_results(); // due to JA

  protected:

    // Data members
    bool                                   my_valid;

    PedigreeDataSet                        my_ped_data;
    model                                  my_model;
    MAXFUN::Results                        my_maxfun_results;
    uint                                   my_subpedigree_count;
    uint                                   my_unconnected_count;
    pg::post_geno_map                      my_type_probs;
    pf::pen_func_map                       my_pfunc_probs;
    
    MlmResidCorrelationCalculator	           my_mlm_resid_corr;
};

//lint -esym(1712,primary_analysis)

/** The primary SEGREG analysis class
 *  This class calculates results of a single SEGREG analysis run, which may
 *  involve one or more models being compared.  The results of this analysis
 *  are returned in the primary_analysis_results object.
 *
 *  Multiple runs may be performed with a single primary analysis object. 
 *  However, the results object will be replaced each time a new analysis is
 *  run.  The references to prior analysis results will therefore become invalid.
 *  Copy any results necessary before the new analysis is performed.
 */
class primary_analysis
{
  public:

    typedef vector<model>                          model_vector;

    primary_analysis(APP::Output_Streams& o);

    ~primary_analysis();

    bool is_quality () const;
    void set_quality(bool);

    const primary_analysis_results& run_analysis
                (const PedigreeDataSet& ped_data,
                 const model_vector&    models,
                 bool                   user_defined = true);
 
    const primary_analysis_results& get_results() const;

    bool check_all_fixed(const PedigreeDataSet& ped_data, const model& some_model);
    vector<std::pair<string,double> > do_single_evaluation(const PedigreeDataSet& ped_data, \
    const model& test_model, double& func_val, vector<string>& par_type);

    static double like_cutoff;

    
  protected:

    /** The evaluation model is a struct for storing the model and its
     *  evaluation tools.  These can then be compared using the comparison
     *  functions.
    */
    class evaluation_model
    {
    public:
    
      evaluation_model(const model& mv, const PedigreeDataSet& ped_data);

      model                mdl;
      MAXFUN::ParameterMgr params;
      segreg_calculator    calculator;
      MAXFUN::Results      current_results;

      bool                is_bad;

    };

    typedef boost::shared_ptr<evaluation_model> eval_model_ptr;
    typedef vector<eval_model_ptr>              eval_model_vector;

    static bool  eval_model_bad (const eval_model_ptr&);
    static bool  model_less_than(const eval_model_ptr&, const eval_model_ptr&);

    class fails_cutoff : public std::unary_function<eval_model_ptr, bool>
    {
    public:
      fails_cutoff(double c) : cutoff(c) { }

      bool operator()(const eval_model_ptr& e) const
      {
        return e->current_results.getFinalFunctionValue() < cutoff; 
      }

    private:

      double cutoff;
    };

    // Unary function which returns true if there are std. errors != 0, false otherwise
    class has_std_err: public std::unary_function<eval_model_ptr, bool>
    {
    public:
      has_std_err() { }

      bool operator()(const eval_model_ptr& e) const
      {
        const MAXFUN::ParameterMgr& params = e->params;
      
        for(int i = 0; i < params.getParamCount(); ++i)
          if(std::abs(params.getParameter(i).getStdError()) > 1e-12) return true;

        return false; 
      }
    };
    
    // Analysis functions
    void                      construct_models(const PedigreeDataSet& ped_data,
                                               const model_vector&    mv);

    double 		      get_initial_penalization(const PedigreeDataSet& ped_data) const;

    double                    get_beta_max() const;

    bool                      test_models(const PedigreeDataSet& ped_data,
                                          bool user_defined_analysis);

    segreg_errors::error_code test_model(uint                   m);

    void                      run_initial_maxfun();
    void                      run_initial_maxfun(uint model_index);
       
    void                      reduce_model_set();

//    void                      run_continuous_final_maxfun();
//    void                      run_discrete_final_maxfun();

    void                      run_final_maxfun();
    void                      run_final_maxfun(uint model_index);

    void                      finalize_model();

    // Output streams
    APP::Output_Streams& out;
    cerrorstream         errors;

    // Data members
    bool my_quality; // 1 = high, 0 = low (faster)

    eval_model_vector        my_models;

    primary_analysis_results my_results;
};

} // End segreg namespace
} // End SAGE namespace

#include "segreg/analysis.ipp"

#endif
