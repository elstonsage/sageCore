#ifndef SEGREG_TRANS_SUB_MODEL_H
#define SEGREG_TRANS_SUB_MODEL_H
//============================================================================
// File:      transm_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/9/01 - created.                              djb
//                                                                          
// Notes:     defines the transmission sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "sub_model_base.h"
#include "freq_sub_model.h"

namespace SAGE
{

namespace SEGREG
{

/// @name transmission sub-model constants
//@{
extern const std::string  TRANSMISSION_NAME;
extern const double       TRANSM_DEFAULT_VALUE;
extern const double       TRANSM_HOMOG_GENERAL_DEFAULT_VALUE;
extern const double       TRANSM_HOMOG_MEND_AA_DEFAULT_VALUE;
extern const double       TRANSM_HOMOG_MEND_AB_DEFAULT_VALUE;
extern const double       TRANSM_HOMOG_MEND_BB_DEFAULT_VALUE;
extern const double       TRANSM_AB_FREE_AA_DEFAULT_VALUE;
extern const double       TRANSM_AB_FREE_BB_DEFAULT_VALUE;
extern const bool         TRANSM_DEFAULT_FIXED;
extern const double       TRANSM_LB;
extern const double       TRANSM_UB;
//@}    

//----------------------------------------------------------------------------
//  Class:    transmission_sub_model
//                                                                          
//----------------------------------------------------------------------------
//lint -esym(1712,transmission_sub_model) <- No default constructor
//
class transmission_sub_model : public SegregSubmodel
{
  friend class model;

  public:
    enum sm_option { no_trans = 1, homog_no_trans, homog_mendelian, homog_general, 
                     general, tau_ab_free, mitochondrial };
                     
    // Constructor/destructor.  
    transmission_sub_model(const genotype_frequency_sub_model* gf_ptr, cerrorstream& errors = sage_cerr);
    transmission_sub_model(const transmission_sub_model& other);
    transmission_sub_model&  operator=(const transmission_sub_model& other);
    virtual ~transmission_sub_model();
    
    // Gets.
    sm_option       option() const;
    virtual string  option_description() const;  
    virtual string  name() const;
    double          prob(genotype_index indiv, genotype_index mother,
                         genotype_index father) const;                        
    bool            default_() const;
    
    bool has_sex_effect() const;
    
    void  dump(std::ostream& out) const;
    
    // Sets.
    bool  set(sm_option opt, const model_input& tau_AA, const model_input& tau_AB, 
              const model_input& tau_BB, bool no_bounds, bool type_missing);
              
    static std::string  option_2_description(sm_option option);
    static std::string  option_2_parameter(sm_option option);
  
    // GCW - We need this public so we can know what they are
    double  tau(genotype_index gt) const;

    /** Tests the submodel for completeness, which is defined as all
     *  taus required are specified. */
    bool is_complete() const;

  protected:

    virtual int  update();
  
    // Option specific set sub-model functions. 
    bool  set_no_trans(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    bool  set_no_trans_hwe(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    bool  set_no_trans_nhwe(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    bool  set_homog_no_trans(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    bool  set_homog_mendelian(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    bool  set_homog_general(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB,
                            double transm_lb, double transm_ub);
    bool  set_general(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB,
                      double transm_lb, double transm_ub);
    bool  set_tau_ab_free(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB,
                          double transm_lb, double transm_ub);
    bool  set_mitochondrial(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    
    // Synchronize w. maxfun.
    int  synchronize_no_trans();
    int  synchronize_homog_no_trans();
    int  synchronize_homog_mendelian();
    int  synchronize_homog_general();
    int  synchronize_general();
    int  synchronize_tau_ab_free();
    int  synchronize_mitochondrial();
    
    // Ancillary functions.
    bool  input_meets_constraints(const model_input& input) const;
    bool  inputs_ok(const model_input& tau_AA, const model_input& tau_AB, const model_input& tau_BB);
    
    void  calc_tau_ab();

    bool  set_homog_general_one_value(const model_input& mi, bool input_is_AA);
    bool  set_homog_general_two_values(const model_input& mi_one,
                                       const model_input& mi_two);
  
    bool  set_general_one_value(const model_input& mi, genotype_index index);
    
    bool  set_tau_ab_free_AB_value(const model_input& tau_AB);
    
    /// Override of the base MAXFUN::Submodel's function which copies
    /// my_taus into my_parameters when the model calls for it.
    ///
    /// \returns 0, always successful.
    virtual int finalizeConfiguration();
    
  private:

    /// Initialize the maxfun initialization matrix.  Sets up three parameters,
    /// with types as defined by the MAXFUN::Parameter::ParamTypeEnum arguments,
    /// and bounds as defined by the lower and upper bound.
    ///
    /// \param AA_type     The initialization type of tau_AA
    /// \param AB_type     The initialization type of tau_AB
    /// \param BB_type     The initialization type of tau_BB
    ///
    /// \param lower_bound The smallest value a transmission can have, generally
    ///                    either TRANSM_LB or NEGATIVE_INF.
    /// \param lower_bound The largest value a transmission can have, generally
    ///                    either TRANSM_UB or POSITIVE_INF.
    void  initialize_parameters(MAXFUN::Parameter::ParamTypeEnum AA_type,
                                MAXFUN::Parameter::ParamTypeEnum AB_type,
                                MAXFUN::Parameter::ParamTypeEnum BB_type,
                                double                           lower_bound,
                                double                           upper_bound);

    /// Copies my_taus into my_parameters as part of finalizing the model
    /// prior to maximization.  Also used at the end of set fuctions for
    /// synchronizing my_parameters with my_taus.
    void finalize_parameters();

    /// Copies the frequency of allele A into my_taus.  Used by the homog_no_trans
    /// model.  Note that the transmission sub model \b must be updated \b after
    /// the frequency sub model for this to work properly.  Always add in
    /// order of least dependent to most dependent sub model.
    void  copy_freq_a_into_taus();
  
    // Data members.  
    sm_option  my_option;
    double     my_taus[NUM_OF_TYPES];
    const genotype_frequency_sub_model* my_gf_ptr;    // For access to frequency of allele A.
    bool  my_default;                                 // True means user has not specified any estimable parameters.
};

} 
} 

#include "segreg/transmission_sub_model.ipp"

#endif
