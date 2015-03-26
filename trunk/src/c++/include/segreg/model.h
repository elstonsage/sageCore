#ifndef SEGREG_MODEL_H
#define SEGREG_MODEL_H
//============================================================================
// File:      model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/14/01 created         -djb
//            8/08/01 Model moved out -gcw
//                                                                          
// Notes:     Defines class, parser, for parsing segreg_analysis block of a
//            SAGE parameter file.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <string>
#include <limits>
#include <bitset>
#include "rped/rped.h"
#include "segreg/sub_model_base.h"
#include "segreg/type_sub_model.h"
#include "segreg/freq_sub_model.h"
#include "segreg/resid_sub_model.h"
#include "segreg/transmission_sub_model.h"
#include "maxfunapi/TransformationSubmodel.h"
#include "segreg/cov_sub_model.h"
#include "segreg/fpmm_sub_model.h"
#include "segreg/onset_sub_model.h"
#include "segreg/ascertainment_sub_model.h"
#include "segreg/prev_sub_model.h"
#include "segreg/definitions.h"
namespace SAGE
{

namespace SEGREG
{

std::string  model_class_2_string(model_class mc);
std::string  primary_type_2_string(primary_type pt);

//----------------------------------------------------------------------------
//  Struct:    model
//                                                                          
//  Purpose:  repository for SEGREG sub-models and other parameter file 
//            information.
//                                                                          
//----------------------------------------------------------------------------
//
enum inconsistency  
{
  i_non_specific,
  i_transm_one_type,            ///< transm. opt. not no_trans or homog_no_trans. and 1 type 
  i_transm_nhwe,                ///< homogeneous opts. in transm. and nhwe. 
  i_prim_cov,                   ///< primary trait also a covariate 
  i_transm_pen_out,             ///< transm is homog_mendelian & pen_func_out is true 
  i_type_prob_one_type,         ///< type_prob. = true & 1 type. 
  i_freq_not_none_one_type,     ///< geno freq sm option not NONE and one type mean 
  i_mv_types,                   ///< mean/var options not compatible in terms of number of types 
  i_freq_nhwe_two_types,        ///< geno freq sm nhwe w two means 
  i_transm_no_trans_two_types,  ///< transm sm no_trans w two means 
  i_binary_nt_no_residuals,     ///< binary trait with no transmission and no residuals invalid
  i_mc_cov_non_exclusive,       ///< The mean and composite trait covariates must not share traits
  i_ms_cov_non_exclusive,       ///< The mean and susc covariates must not share traits in age-of-onset models
  i_mv_cov_non_exclusive,       ///< The mean and variance covariates traits must not share traits in age-of-onset models
  i_sc_cov_non_exclusive,       ///< The susc and composite trait coviarates must not share traits in age-of-onset models
  i_num_i                       ///< number of inconsistencies 
};

class model
{
  public:
    friend class parser;

    model(cerrorstream& errors = sage_cerr);
    model(const model& other);
    model&  operator=(const model& other);
    void reset(cerrorstream& errors);
    
    std::bitset<i_num_i>  check_consistency();
    
    // Sub-models -- Ordered to reflect dependencies
    genotype_specific_mean_sub_model            mean_sub_model;
    genotype_specific_variance_sub_model        var_sub_model;
    genotype_specific_susceptibility_sub_model  susc_sub_model;

    CompositeTraitSubmodel                      comp_trait_sub_model;

    MeanCovariateSubmodel                       mean_cov_sub_model;
    VarianceCovariateSubmodel                   var_cov_sub_model;
    SusceptibilityCovariateSubmodel             susc_cov_sub_model;

    finite_polygenic_mixed_model_sub_model      fpmm_sub_model;
    onset_sub_model                             ons_sub_model;
    residual_correlation_sub_model              resid_sub_model;
    MAXFUN::TransformationSubmodel              transf_sub_model;
    genotype_frequency_sub_model                freq_sub_model;
    transmission_sub_model                      transm_sub_model;
    ascertainment_sub_model                     ascer_sub_model;
  
    prevalence_sub_model                        prev_sub_model;
    
    /// Stores options regarding maxfun's debug features
    ///
    MAXFUN::DebugCfg                            my_maxfun_debug;

    // Model qualifiers
    
    /** the type_dependent_sub_model returns either the mean or susc submodel,
     *  dependent upon the primary trait type as follows:
     *
     *   - pt_CONTINUOUS -> mean_sub_model
     *   - pt_BINARY     -> susc_sub_model
     *   - pt_ONSET      -> dependent on which is type dependent.  See onset_sub_model
     */
    const genotype_specific_mean_susc_sub_model& type_dependent_sub_model() const;
          genotype_specific_mean_susc_sub_model& type_dependent_sub_model();

    // Accessors.
    std::string  get_file_name_root()   const;
    std::string  get_title()            const;
    
    model_class     get_model_class()         const;
    primary_type    get_primary_trait_type()  const;
    bool            get_each_pedigree()       const;
    bool            get_pen_func_output()     const;
    bool            get_type_prob()           const;
    std::string     get_primary_trait()       const;
    bool            get_type_missing()        const;
    bool            get_trans_missing()       const;

    void         set_primary_trait(const std::string& trait_name)
                 { primary_trait = trait_name; }
    void         set_primary_trait_type(primary_type p)
                 { primary_trait_type = p; }
    
    bool has_sex_effect() const;
    
    static set<FPED::MemberConstPointer> cond_mem_set; // for members of conditioned subset due to JA 
  private:
  
    /// Makes certain the covariate sub models of the model point to the
    /// right things.
    void rebuild_covariates();

    // - Meta-constraints.
    //
    bool  transm_one_type()           const;
    bool  transm_nhwe()               const;
    bool  prim_cov()                  const;
    bool  transm_pen_out()            const;
    bool  type_prob_one_type()        const;
    bool  freq_not_none_one_type()    const;
    bool  mv_types()                  const;
    bool  mv_bad_mean_values();                 // Uses sub_model_base::get_parameters wh. is non-const.
    bool  mv_bad_mean_bounds();
    bool  mv_bad_mean_statuses();
    bool  freq_nhwe_two_types()       const;
    bool  transm_no_trans_two_types() const;
    bool  binary_nt_no_residuals()    const;

    void  test_covariate_exclusivity(bitset<i_num_i>& incon) const;
    
    // - Meta-constraint fixes.
    //
    void  fix_freq_not_none_one_type();
    void  fix_mv_bad_mean_bounds();
    void  fix_mv_bad_mean_statuses();
    void  fix_mv_inconsistencies();
  
    // Data members.
    std::string  title;
    std::string  file_name_root;
    
    // Misc. user options.
    model_class     m_class;
    bool            each_pedigree;
    bool            pen_func_output;
    bool            type_prob;
    
    // - User asks for a commingling analysis by omitting a mean sub-model.
    //   Similarly, user asks for a transmission analysis by specifying a mean
    //   sub-model and not a transmission sub-model.
    //
    bool            mean_missing;   //< Mean sub model was missing when parsed
    bool            susc_missing;   //< Susc sub model was missing when parsed
    bool            trans_missing;
    
    // Other parameters.
    std::string    primary_trait;
    primary_type   primary_trait_type;
};

// bool SEGREG::model::no_poly = false;
}
}

#include "segreg/model.ipp"

#endif
