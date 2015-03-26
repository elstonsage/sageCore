#ifndef ITERATIVE_MODELS_H
#define ITERATIVE_MODELS_H

#include "app/output_streams.h"
#include "segreg/mean_split.h"
#include "segreg/type_sub_model.h"
#include "segreg/model.h"
#include "segreg/analysis.h"

namespace SAGE
{
namespace SEGREG
{

/// When maximizing models in the absence of initial estimates, SEGREG uses an
/// iterative process to develop the model from simpler models.  In most
/// cases this works quite well.  In this process, we pass through a matrix
/// of models using different mean and transmission options.  Each is stored
/// for later retrieval along with maxfun results from its maximization.
///
/// The process here is quite simple.  When a target_model is given, all
/// elements are cleared.  Then, as requests are given for certain other
/// models, they are compared with the target.  Where features are missing
/// in the target or are different from the requested model, then the model
/// goes to models which are heirarchically lower than it for its initial
/// estimates, which in turn recursively calls models lower than it as
/// needed, eventually leading to the lowest level model that it has not
/// previously computed.  The chain is then unrolled, building the final
/// model from multiple levels of heirarchy.
///
/// The only real restriction to this process is that variance models are
/// restricted to the 'one' option.  Variances can be computed externally
/// from the results of this complex generation of models.

class Iterative_Models
{
  public:
    
    typedef genotype_specific_mean_sub_model::sm_option mean_option;
    typedef transmission_sub_model::sm_option           transm_option;
    typedef genotype_frequency_sub_model::sm_option     freq_option;

    Iterative_Models(const PedigreeDataSet& ped_data,
                     APP::Output_Streams& output,
                     primary_analysis& a);

    void         set_target_model(const model& target);
    const model& get_target_model() const;

    void  set_output(bool);
    bool  get_output() const;

    void  set_output_mean(bool);
    bool  get_output_mean() const;

    void  set_output_transm(bool);
    bool  get_output_transm() const;
    
    /// The main function to this class, returns the appropriate model or an
    /// empty results structure if it cannot be generated.  It is
    /// intelligent, modifying the mean, transmission combination where they
    /// are non-sensical (ie, mean = two, transm = no_trans returns the two
    /// mean, transm = homog_no_trans, as they are equivalent)

    const primary_analysis_results& get_model(mean_option, transm_option,
                                              bool quality = true);
  
    void set_one_mean_res(double ); // (due to JA for new & improved initial estimates)
    void set_one_susc_res(double ); // (due to JA for new & improved initial estimates)

    static vector<primary_analysis_results> intermax;
    // results of intermediate maximization
  protected:

    typedef primary_analysis::model_vector model_vector;

    struct model_results;

    /// Clear clears the map and makes all the models into a base model.
    void clear();

    /// Simplify the model for processing, if necessary.
    void create_initial_model(model&);

    bool create_initial_type   (model&);
    bool create_initial_var    (model&, double variance = QNAN);
    bool create_initial_transm (model&);
    bool create_initial_freq   (model&);
    bool create_initial_fpmm   (model&);
    bool create_initial_onset  (model&);

    /// Get the initial (simplest possible) model
    const model& get_initial_model();

    /// Boolean comparisons versus the user specified target
    //@{

    /// Returns true if the target mean model is complete and the mean option == it
    bool use_target_mean_model(mean_option mean_opt) const;

    /// Returns true if the target freq model is complete and is hwe if
    /// mean_opt is not a three model.
    bool use_target_freq_model(mean_option mean_opt) const;

    /// Returns true if the target transm model is complete and the transm option == it
    bool use_target_transm_model(transm_option transm_opt) const;

    //@}

    // Various update routines

    /// update all models in v to the target_mean_model
    void update_means_to_target(model_vector& v) const;

    /// update all models in v to the target_freq_model
    void update_freqs_to_target(model_vector& v) const;

    /// transform the transmission model of elements of v to homog_general,
    /// and expand vector as necessary
    void transform_models_to_homog_general(model_vector& v) const;

    /// transform the transmission model of elements of v to tau_ab_free,
    /// and expand as necessary
    void transform_models_to_tau_ab_free  (model_vector& v) const;

    /// transform the transmission model of elements of v to general
    void transform_models_to_general      (model_vector& v) const;

    void set_two_variance(model& m);

    void set_three_variance(model& m, bool use_freq);

    void create_two_mean_estimates  (model_vector& v, mean_option m);
    void create_three_mean_estimates(model_vector& v, mean_option m) const;

    void create_continuous_two_mean_estimates  (model_vector& v, mean_option m);
    void create_discrete_two_mean_estimates    (model_vector& v, mean_option m);

    void get_three_common_models  (model_vector& v, mean_option m, transm_option t);

    void transform_three_means(model_vector& v, mean_option m) const;

    /// Calculation Routines
    //@{
    
    template<int mean_count, transm_option T>
    const primary_analysis_results& get_model_results(mean_option m, bool quality);

    template<int mean_count, transm_option T>
    void get_starting_models(model_vector& v, mean_option m);

    template<int mean_count, transm_option T>
    void transform_starting_models(model_vector& v, mean_option m);

    //@}

    /// Switchboard operations
    //@{
    
    const primary_analysis_results& two  (mean_option m, transm_option t, bool q);
    const primary_analysis_results& three(mean_option m, transm_option t, bool q);

    //@}

    /// Perform analysis

    void run_models(mean_option m, transm_option t, model_results& ret, const model_vector& v, bool q);

    void produce_model_failure(mean_option m, transm_option t, model_results& res);


    // Data Members

    model my_target_model;
    model my_initial_model;

    bool  initial_model_built;    //< Has the initial model been built?
    bool  target_mean_complete;   //< Is the target mean model complete?
    bool  target_transm_complete; //< Is the target transm model complete?
    bool  target_freq_complete;   //< Is the target freq model complete?

    struct model_results
    {
      /// The availaibility_t type classifies the state of the results
      
      enum availability_t
      {
        empty,       //< No results yet.  This is the initial state.
        good,        //< Results are good!  Can be used as a starting state in later calculations.
        bad          //< Results are bad.  Should not be used as a starting state in later calculations.
      };
    
      model_results();
          
      void clear();                //< Clears the results and sets to empty
              
      bool is_good       () const; //< returns true if good
      bool is_bad        () const; //< returns true if bad
      bool is_empty      () const; //< returns true if empty
      bool is_available  () const; //< returns true if !empty

      availability_t           availability;
      primary_analysis_results result;
    };
    
    /// Models are classified by their type and transmission options.
    /// Thus we use a simple class to compare them.
    
    struct model_classification
    {
      model_classification();
      model_classification(mean_option m, transm_option t);

      mean_option   mean_opt;
      transm_option transm_opt;
      
      bool operator< (const model_classification& rhs) const;
    };

    std::map<model_classification, model_results> my_iterative_models;

    mean_split my_trait_sample;

    // non data storage elements

    //lint -e1725
    PedigreeDataSet        my_ped_data;
    primary_analysis&      analyzer;
    APP::Output_Streams&   out;
    //lint +e1725

    bool output;
    bool output_mean;
    bool output_transm;
    
    double one_mean_res; // (due to JA for new & improved initial estimates)
    double one_susc_res; // (due to JA for new & improved initial estimates)

    void output_max_line(mean_option m, transm_option t, bool done = false) const;
};

}
}

#include "segreg/iterative_models.ipp"

#endif
