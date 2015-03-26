#ifndef SEGREG_RESID_SUB_MODEL_H
#define SEGREG_RESID_SUB_MODEL_H
//============================================================================
// File:      resid_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/2/01 - created.                                   djb
//                                                                          
// Notes:     defines the residual correlations sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "error/internal_error.h"
#include "sub_model_base.h"
#include "definitions.h"

namespace SAGE
{

namespace SEGREG
{

/// @name resid sub-model constants
//@{
const int                 NUM_OF_CORRS = 8;

extern const std::string  RESID_NAME;
extern const double       RESID_SP_DEFAULT_VALUE;
extern const double       RESID_DEFAULT_VALUE;
extern const double       RESID_EPSILON;
extern const bool         RESID_DEFAULT_FIXED;           
extern const double       RESID_LB;                       ///< Per rce 8/14/01.
extern const double       RESID_UB;
//@}

//----------------------------------------------------------------------------
//  Class:    residual_correlation_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class residual_correlation_sub_model : public SegregSubmodel
{
  public:
    enum sm_option { equal_po_ss, equal_po, arb };
    enum corr  { fm, fs, ms, fd, md, bb, bs, ss };

    /// The indices define the index into my_parameters based upon the
    /// sm_option.  They are prefaced by a code indicating the sm_option:
    /// eps = equal_po_ss, ep = equal_po, arb = arb.
    
    enum indices { eps_fm = 0, eps_po_ss,
                   ep_fm  = 0, ep_po,     ep_ss,
                   arb_fm = 0, arb_fo,    arb_mo, arb_ss };

    // Constructor/destructor.  
    residual_correlation_sub_model(cerrorstream& errors = sage_cerr);
    residual_correlation_sub_model(const residual_correlation_sub_model& other);
    residual_correlation_sub_model&  operator=(const residual_correlation_sub_model& other);
    virtual ~residual_correlation_sub_model();
    
    // Gets.
    sm_option       option() const;
    virtual string  option_description() const;  
    virtual string  name() const;

    bool            is_in_default_mode() const;
    
    double          father_mother_correlation() const;
    double          mother_son_correlation() const;
    double          father_son_correlation() const;
    double          mother_daughter_correlation() const;
    double          father_daughter_correlation() const;
    double          brother_brother_correlation() const;
    double          sister_sister_correlation() const;
    double          brother_sister_correlation() const;
    
    double          correlation(corr rs) const;
    bool            is_correlation_fixed(corr rs) const;

    /// Returns true if there are residuals which can be != 0.
    bool            has_residuals() const;
    
    /// Returns \c true if there 
    bool            has_sex_effect() const;
    
    double          alpha_mother(double tf, double tm) const;   //< Equ. 52.
    double          alpha_father(double tf, double tm) const;   //< Equ. 53.
    double          delta       (double tf, double tm) const;   //< Equ. 55.
    
    double          alpha_mother(bool tf, bool tm) const;   //< Equ. 52.
    double          alpha_father(bool tf, bool tm) const;   //< Equ. 53.
    double          delta       (bool tf, bool tm) const;   //< Equ. 55.
    
    void  dump(std::ostream& out) const;
    
    // Sets.
    bool  set_as_default(model_class m_class);
    
    bool  set(sm_option opt, model_input fm, const model_input mo, 
                             model_input fo, const model_input ss,
                             model_class m_class);                         
//    bool  set(sm_option opt, const model_input& fm, const model_input& ms, 
//                             const model_input& fs, const model_input& md,
//                             const model_input& fd, const model_input& bb,
//                             const model_input& ss, const model_input& bs );   // NOT implemented.
              
//    void  supply_missing_values(model_class m_class);

    static std::string  option_2_description(sm_option option);
    static std::string  option_2_parameter(sm_option option);
  
    void  calculate_alpha_and_delta();

  protected:

    void  supply_missing_values();

    virtual int update();
    
    // Option specific synchronization.
    bool  synchronize_equal_po_ss();
    bool  synchronize_equal_po();
    bool  synchronize_arb();
    
    // Internal synchronization.
    void  internally_synchronize();
    void  internally_synchronize_equal_po_ss();
    void  internally_synchronize_equal_po();
    void  internally_synchronize_arb();
    
    /// Calculate the alpha and deltas, which is typically what we really
    /// want


    // Option specific sets.
    bool  set_equal_po_ss(const model_input& fm, const model_input& mo, const model_input& fo, const model_input& sib_sib);
    bool  set_equal_po(const model_input& fm, const model_input& mo, const model_input& fo, const model_input& sib_sib);
    bool  set_arb(const model_input& fm, const model_input& mo, const model_input& fo, const model_input& sib_sib);

    virtual int finalizeConfiguration();
    
    /// The initialize function is called after all the initial correlations
    /// are set.  It sets all the other (dependent) correlations and builds
    /// the maxfun parameters using the other initialize_* functions as helpers.
    void  initialize();
    
    /// Initializes all the dependent correlaitons for the equal_po_ss model and
    /// builds the maxfun parameters
    void  initialize_equal_po_ss();

    /// Initializes all the dependent correlaitons for the equal_po model and
    /// builds the maxfun parameters
    void  initialize_equal_po();

    /// Initializes all the dependent correlaitons for the arb model and
    /// builds the maxfun parameters
    void  initialize_arb();
    
    // Ancillary functions.
    bool  input_meets_constraints(model_input& input, double epsilon = 0);
    bool  set_po_ss(const model_input& mo, const model_input& fo, const model_input& sib_sib);

    void  supply_missing_values(double pov, double ssv);

    double get_lower_bound() const;
    double get_upper_bound() const;
    
    // Data members.  
    sm_option   my_option;
    
    model_class my_model_class;

    double     my_corrs[NUM_OF_CORRS];
    bool       my_corrs_fixed[NUM_OF_CORRS];

    bool       my_is_in_default_mode;  ///< True means user has not specified a residual model

    // Storage of equation 52 - 55 values.
    
    double my_alpha_mother[4];
    double my_alpha_father[4];
    double my_delta[4];
};

} 
} 

#include "segreg/resid_sub_model.ipp"

#endif
