#ifndef SEGREG_PREV_SUB_MODEL_H
#define SEGREG_PREV_SUB_MODEL_H
//============================================================================
// File:      prev_sub_model.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   2004-03-18  created.                              gcw
//                                                                          
// Notes:     Defines and impelents prevalence constraints and estimates for SEGREG
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "numerics/log_double.h"
#include "sub_model_base.h"
#include "type_sub_model.h"
#include "cov_sub_model.h"
#include "freq_sub_model.h"
#include "fpmm_sub_model.h"
#include "onset_sub_model.h"

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace SAGE
{
namespace SEGREG
{

/// \name Prevalence Calculation Objects
///
/// The prevalence sub model stores data for two sets of calculations,
/// prevalence estimates and prevalence constraints.  These are very similar,
/// but there are important differences.
///
/// A prevalence estimate is an estimate of the prevalence of a binary (or
/// onset) trait given a set of covariates drawn from the susceptibililty
/// covariates (or mean or variance covariates in the case of onset traits)
/// with associated values. Thus the data storage is primarily that of the
/// covariates and their values.  From this data, we calculate an estimated
/// prevalence for each maximization step (creating a dependent parameter in
/// MAXFUN).
///
/// Prevalence constraints take this calculation one step further,
/// calculating a likelihood penalty based on the estimated prevalence.  To
/// do this, it requires two numbers, the sample size (N) and the number
/// affected within that sample (R).  The resulting penalty directly modifies
/// the MAXFUN likelihood.  It does not, however, create a variable within
/// MAXFUN.
///
/// In addition, in either case, if the trait is age of onset, we require an
/// age at which the estimate/constraint is to be calculated.
///
/// The prevalence_sub_model is primarily two vectors, one of estimates and
/// one of constraints.  Each element in these vectors must be unique, and
/// must also have valid data.  This makes building the prevalence_sub_model
/// 'in place' difficult.
///
/// Therefore, the prevalence_sub_model is not built with set() functions as
/// other sub_models, but rather with a "builder" object.  The builder
/// provide an iterative method of creating a new prevalence
/// estimate/constraint to be added to the sub model.  This allows much
/// smoother construction, as errors can be reported as each phase is
/// complete.  As a final step, once an estimate/constraint is finished, it
/// can then be added to the sub model with the add_estimate() or
/// add_constraint() function.  The associated builder is then cleared, so
/// that it can be reused for another constraint/estimate, or discarded.
//@{

/// The psm_info class stores all information which is necessarily common to
/// the prevalence_sub_model and the builder objects.
class psm_info
{
protected:

    friend class psm_builder;
    friend class prevalence_sub_model;

    /// Defines a prevalence covariate.  May be used for either prevalence
    /// constraints and/or prevalence estimates

    struct prev_cov
    {
      prev_cov();
      
      bool operator==(const prev_cov&) const;
                            
      std::string  my_name;
      double       my_value;

      size_t       my_index;
    };

    typedef std::vector<prev_cov> prev_cov_vector;

    /// Defines a prevalence estimate. This involves a set of prev_covariate
    /// items and the age (for onset traits).
    
    struct prev_estimate
    {
      bool operator==(const prev_estimate&) const;

      /// Susceptibility covariates
      prev_cov_vector my_susc_covs;

      /// Mean Covariates (only for age of onset)
      prev_cov_vector my_mean_covs;

      /// Variance Covariates (only for age of onset)
      prev_cov_vector my_var_covs;

      double          my_number_affected; ///< R
      double          my_sample_size;     ///< N

      double          my_age;             ///< Age (only for age of onset)
    };

    /// Defines a prevalence constraint.  This is identical to the
    /// prevalence estimate with two additional parameters: R (number
    /// affected) and N (sample size) values for population prevalence.

    struct prev_constraint : public prev_estimate
    {
      bool operator==(const prev_constraint&) const;

      double          my_number_affected; ///< R
      double          my_sample_size;     ///< N
    };

};

/// The psm_builder is used to construct and store the information which will
/// build a prevalence estimate or constraint.
class psm_builder
{
public:

    friend class prevalence_sub_model;

    psm_builder(const MeanCovariateSubmodel*            mean_covs,
                const SusceptibilityCovariateSubmodel*  susc_covs,
                const VarianceCovariateSubmodel*        var_covs);
            
    /// When covariates are added, the add_covariate function returns one of
    /// the following states, which must be evaluated by the application.
    enum add_cov_ret_type
    {
      valid,     ///<  Valid.
      duplicate, ///<  Covariate already in set, and replacement done.  Assumes valid.
      unknown,   ///<  Covariate name not found.
      value_err  ///<  Value is not valid (not finite).
    };

    /// Adds a covariate to the list with a value for that covariate. 
    /// Covariate must be present in the mean covariates, present in the
    /// susceptibility covariates, or present in the variance covariates. 
    /// Covariates may be added multiple times, but values override prior
    /// values. 
    ///
    /// Returns one of the add_cov_ret_type above.
    
    add_cov_ret_type add_covariate(const string& cov_name, double value);

    /// Sets the age.  Returns 0 if valid, 1 if the value is NaN or is <= 0.

    size_t set_age(double age);

    /// Sets the sample size.  Returns 0 if valid, 1 if the value is not
    /// finite or is <= 0.  Does not cross check with number affected, which
    /// may not be set yet.  Also does not care, at this point, if the
    /// object being built is a constraint or an estimate (the latter does
    /// not require sample size or number affected).

    size_t set_sample_size     (double sample_size);

    /// Sets the number_affected.  Returns 0 if valid, 1 if the value is not
    /// finite or is <= 0.  Does not cross check with sample size, which may
    /// not be set yet.  Also does not care, at this point, if the object
    /// being built is a constraint or an estimate (the latter does not
    /// require sample size or number affected).

    size_t set_number_affected (double number_affected);

protected:

   /// Clear the prevalence data.
   
   void clear_prev_data();

   /// add_covariate() helper functions
   //@{

   add_cov_ret_type modify_value (double& lhs, double value);
   
   add_cov_ret_type add_mean_cov (const string& cov_name, double value);
   add_cov_ret_type add_susc_cov (const string& cov_name, double value);
   add_cov_ret_type add_var_cov  (const string& cov_name, double value);

   //@}
  

   /// Data for verification of added covariates.
   //@{

   const MeanCovariateSubmodel*            my_mean_covs;
   const SusceptibilityCovariateSubmodel*  my_susc_covs;
   const VarianceCovariateSubmodel*        my_var_covs;
   //@}

   /// The data for the estimate/constraint.
   
   /// We actually store values for all means and susceptibilities.  This
   /// is to make it easier and more efficient to build these objects
   /// consistently.  When stored in the prevalence_sub_model, it is done
   /// more efficiently.
   ///
   /// These are the elements cleared by the clear_prev_data() function.
   //@{
   
   vector<double> my_mean_cov_values;
   vector<double> my_susc_cov_values;
   vector<double> my_var_cov_values;
   
   double my_age;

   double my_sample_size;
   double my_number_affected;

   //@}

};

/// The prevalence_sub_model is what is actually used by the SEGREG model.  It
/// stores the set of prevalence estimates and prevalence constraints.
class prevalence_sub_model : public SegregSubmodel
{
  public:

    friend class model;

    prevalence_sub_model(const genotype_specific_mean_sub_model*           means,
                         const genotype_specific_susceptibility_sub_model* suscs,
                         const genotype_specific_variance_sub_model*       vars,
                         const MeanCovariateSubmodel*                   mean_covs,
                         const SusceptibilityCovariateSubmodel*         susc_covs,
                         const VarianceCovariateSubmodel*               var_covs,
                         const genotype_frequency_sub_model*               freqs,
                         const finite_polygenic_mixed_model_sub_model*     fpmms,
                         const onset_sub_model*                            onsets,
                         const MAXFUN::TransformationSubmodel*             transforms,
                         bool                                              fpmm_option,
                         bool                                              onset_option);

    prevalence_sub_model(const prevalence_sub_model&);

    prevalence_sub_model& operator= (const prevalence_sub_model& psm);

    ~prevalence_sub_model();

    virtual string name() const;
    virtual string option_description() const;
    
    /// Modify the fpmm option.  Return true if the fpmm option was set,
    /// false otherwise.  The only reason for this to fail is if the onset
    /// option is currently true (which requires fpmm)
    
    bool set_fpmm_option(bool using_fpmm);

    /// Modify the onset option.  Return true if the onset option was set,
    /// false otherwise.  Note that the only reason for the onset option to
    /// fail would be that there are constraints and/or estimates already
    /// existing in the current model or that the fpmm option is currently
    /// false.
    
    bool set_onset_option(bool using_age_onset);

    // Build operations

    /// Return types when elements are added. 
    
    enum add_elt_ret_type
    {
      valid,         ///< The standard valid value
      na_not_finite, ///< Number Affected not finite
      ss_not_finite, ///< Sample Size not finite
      na_gt_ss,      ///< Number Affected >= Sample Size
      age_required,  ///< Age required but no age specified
      age_present,   ///< Age present when it's not needed
      duplicate      ///< Element already present
    };
   
    /// Adds a constraint.  If valid, the psm_builder object is cleared and
    /// 0 returned.  Returns error state if there is a problem.
     
    add_elt_ret_type add_constraint(psm_builder&);

    /// Adds a estimate.  If valid, the psm_builder object is cleared and 0
    /// returned.  Returns error state if there is a problem.
     
    add_elt_ret_type add_estimate(psm_builder&);

    // Get options
    
    /// Returns the log_double product of all the prevalence constraint
    /// penalties.
    log_double get_prevalence_penalty() const;

    /// Returns number of prevalence estimates
    size_t get_estimate_count() const;

    /// Returns the age for onset estimates
    double get_estimate_age (size_t est) const;

    /// Returns number of susceptibility covariates in a particular estimate
    size_t get_estimate_susc_covariate_count(size_t) const;

    /// Returns the name of a susceptibility covariate in an estimate
    const string& get_estimate_susc_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a susceptibility covariate in an estimate
    double get_estimate_susc_covariate_value(size_t est, size_t cov) const;

    /// Returns number of mean covariates in a particular estimate
    size_t get_estimate_mean_covariate_count(size_t) const;

    /// Returns the name of a mean covariate in an estimate
    const string& get_estimate_mean_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a mean covariate in an estimate
    double get_estimate_mean_covariate_value(size_t est, size_t cov) const;

    /// Returns number of variance covariates in a particular estimate
    size_t get_estimate_var_covariate_count(size_t) const;

    /// Returns the name of a variance covariate in an estimate
    const string& get_estimate_var_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a variance covariate in an estimate
    double get_estimate_var_covariate_value(size_t est, size_t cov) const;

    /// Returns number of prevalence constraints
    size_t get_constraint_count() const;

    /// Returns the age for onset constraints
    double get_constraint_age (size_t) const;

    double get_constraint_number_affected(size_t) const;
    double get_constraint_sample_size    (size_t) const;

    /// Returns number of susceptibility covariates in a particular constraint
    size_t get_constraint_susc_covariate_count(size_t) const;

    /// Returns the name of a susceptibility covariate in an constraint
    const string& get_constraint_susc_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a susceptibility covariate in an constraint
    double get_constraint_susc_covariate_value(size_t est, size_t cov) const;

    /// Returns number of mean covariates in a particular constraint
    size_t get_constraint_mean_covariate_count(size_t) const;

    /// Returns the name of a mean covariate in an constraint
    const string& get_constraint_mean_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a mean covariate in an constraint
    double get_constraint_mean_covariate_value(size_t est, size_t cov) const;

    /// Returns number of variance covariates in a particular constraint
    size_t get_constraint_var_covariate_count(size_t) const;

    /// Returns the name of a variance covariate in an constraint
    const string& get_constraint_var_covariate_name (size_t est, size_t cov) const;

    /// Returns the value used by a variance covariate in an constraint
    double get_constraint_var_covariate_value(size_t est, size_t cov) const;

protected:

    /// build helpers

    add_elt_ret_type add_constraint_internal (const psm_builder&);
    add_elt_ret_type add_estimate_internal   (const psm_builder&);
    
    /// Moves prevalence covariates from the builder into the constraint/estimate 
    void copy_prev_covs(psm_info::prev_estimate&, const psm_builder&);

    virtual int  update();

    // Calculation routines

    /// Calculates adjustment to base susc due to covariates
    double calculate_susc_adj(const psm_info::prev_cov_vector& susc_covs) const;

    /// Calculates adjustment to base mean due to covariates
    double calculate_mean_adj(const psm_info::prev_cov_vector& mean_covs) const;

    /// Calculates adjustment to base var due to covariates
    double calculate_var_adj(const psm_info::prev_cov_vector& var_covs) const;

    /// Calculates an adjustment due to covaraites from the sub_model given
    double calculate_adj(const psm_info::prev_cov_vector& covs,
                         const CovariateSubmodel*     csm) const;

    /// Calculates susceptibility
    double calculate_susc(genotype_index u, const psm_info::prev_cov_vector& susc_covs) const;

    /// Calculates mean
    double calculate_mean(genotype_index u, const psm_info::prev_cov_vector& mean_covs) const;

    /// Calculates var
    double calculate_var (genotype_index u, const psm_info::prev_cov_vector& var_covs ) const;

    /// Calculates alpha: \f$\alpha = frac{ \pi } { \sqrt{3} \sigma }
    double calculate_alpha(genotype_index u, const psm_info::prev_cov_vector& var_covs) const;

    /// Calculates transformed age
    double calculate_age(const psm_info::prev_estimate&) const;

    /// Calculates the penetrance given main type u
    double penetrance(const psm_info::prev_estimate&, genotype_index u) const;

    /// Calculates the penetrance given main type u and polygenic value v
    double penetrance(const psm_info::prev_estimate&, genotype_index u, size_t v) const;

    /// Calculates P_rev, the sum of the calculated penetrances * the
    /// frequencies
    double calculate_prev           (const psm_info::prev_estimate&) const;

    /// Calculates the prevalence penalty to the likelihood.
    log_double calculate_prev_penalty(const psm_info::prev_constraint&) const;

    /// Models needed for this to work
    //@{

    const genotype_specific_mean_sub_model*                my_means;
    const genotype_specific_susceptibility_sub_model*      my_suscs;
    const genotype_specific_variance_sub_model*            my_vars;
    const MeanCovariateSubmodel*                        my_mean_covs;
    const SusceptibilityCovariateSubmodel*              my_susc_covs;
    const VarianceCovariateSubmodel*                    my_var_covs;
    const genotype_frequency_sub_model*                    my_freqs;
    const finite_polygenic_mixed_model_sub_model*          my_fpmms;
    const onset_sub_model*                                 my_onsets;
    const MAXFUN::TransformationSubmodel*                  my_transforms;

    /// True if fpmm is used.
    bool my_fpmm_option;

    /// True if age of onset is used.
    bool my_onset_option;

    //@}

    typedef std::vector<psm_info::prev_estimate>   prev_estimate_vector;
    typedef std::vector<psm_info::prev_constraint> prev_constraint_vector;
            
    /// Data for the prevalence estimates
    prev_estimate_vector   my_estimate_data;

    /// Data for the prevalence constraints
    prev_constraint_vector my_constraint_data;
};

//@}
}
}

#include "segreg/prev_sub_model.ipp"

#endif
