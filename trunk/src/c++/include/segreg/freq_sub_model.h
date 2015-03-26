#ifndef SEGREG_FREQ_SUB_MODEL_H
#define SEGREG_FREQ_SUB_MODEL_H
//============================================================================
// File:      freq_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/9/01  - created.                          djb
//                                                                          
// Notes:     defines the genotype frequency sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "segreg/sub_model_base.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/type_sub_model.h"
#include "error/internal_error.h"

namespace SAGE
{

namespace SEGREG
{
/// @name genotype frequency sub-model constants
//@{
extern const std::string  GENO_FREQ_NAME;
extern const double       PROB_AA_DEFAULT_VALUE;
extern const double       PROB_AB_DEFAULT_VALUE;
extern const double       PROB_BB_DEFAULT_VALUE;
extern const double       PROB_AA_DEFAULT_VALUE_NONE;         ///< Used w. 'none' option.
extern const double       PROB_AB_DEFAULT_VALUE_NONE;         ///< Used w. 'none' option.
extern const double       PROB_BB_DEFAULT_VALUE_NONE;         ///< Used w. 'none' option.
extern const double       PROB_TOTAL;
extern const bool         PROBS_FIXED_DEFAULT;
extern const double       FREQ_EPSILON;                       ///< Per rce 7/15/01.  Used to check prob. total.
extern const double       FREQ_A_DEFAULT_VALUE;
extern const double       FREQ_A_DEFAULT_VALUE_NONE;          ///< Used w. 'none' option.
extern const double       FREQ_A_LB;        
extern const double       FREQ_A_UB;
extern const double       FREQ_A_EPSILON;      
extern const double       PROB_LB;        
extern const double       PROB_UB;        
extern const double       CORR_DEFAULT_VALUE;
extern const bool         CORR_DEFAULT_FIXED;           
extern const double       CORR_LB;
extern const double       CORR_UB;
//@}

//----------------------------------------------------------------------------
//  Class:    genotype_frequency_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class genotype_frequency_sub_model : public SegregSubmodel
{
  public:
    static const int NUM_OF_PROBS = 3;
    enum sm_option { hwe = 1, nhwe, NONE };  // 'none' clashes w. enum value in SEGREG namespace.
    
    // Constructor/destructor.  
    genotype_frequency_sub_model(cerrorstream& errors = sage_cerr);
    genotype_frequency_sub_model(const genotype_frequency_sub_model& other);
    genotype_frequency_sub_model&  operator=(const genotype_frequency_sub_model& other);
    virtual ~genotype_frequency_sub_model();
    
    // Gets.
    sm_option       option() const;
    virtual string  option_description() const;  
    virtual string  name() const;
    double          freq_A() const;
    double          prob(genotype_index i) const;

    double          prob(genotype_index i, genotype_index s) const;  // indiv. genotype / spouse genotype   
    bool            default_() const;

    void  dump(std::ostream& out) const;
    
    bool  set(sm_option opt, double freq_A, double prob_AA, 
              double prob_AB, double prob_BB, const model_input& corr, 
              const genotype_specific_mean_susc_sub_model& type_sm, 
              bool type_missing, bool trans_missing, bool probs_fixed = false);
    bool  set(sm_option opt, double freq_A, double prob_AA, 
              double prob_AB, double prob_BB, const model_input& corr, bool probs_fixed = false);
    bool  set_none();
              
    static std::string  option_2_description(sm_option option);
    static std::string  option_2_parameter(sm_option option);
  
    /** Tests the submodel for completeness, which is defined as all
     *  appropriate frequencies being finite. */
    bool is_complete() const;

  protected:

    virtual int update();
    
    // Ancillary set functions.
    bool  set_hwe(double freq_A, double prob_AA, 
                  double prob_AB, double prob_BB, const model_input& corr, bool probs_fixed);
    bool  set_nhwe(double freq_A, double prob_AA, 
                   double prob_AB, double prob_BB, const model_input& corr, bool probs_fixed);
    
                   
    // Option specific synchronization functions (sync w. Maxfun).
    int  synchronize_hwe();
    int  synchronize_nhwe();
    
    // Ancillary functions.
    bool  freq_A_input_meets_constraints(double& input, bool fixed);

    bool  set_one_prob(bool probs_fixed);
    bool  set_two_probs(double prob_one, genotype_index index_one,
                        double prob_two, genotype_index index_two, bool probs_fixed);

    /// \brief my_parameters initialization
    //@{
    /// Creates the my_parameters vector by calling the right initialize
    /// routine for the option.  Always successful, so always returns 0.
    virtual int finalizeConfiguration();
                        
    /// Initializes my_parameters for the none option (0 parameters)
    ///
    void  initialize_none ();

    /// Initializes my_parameters for the hwe option
    ///
    void  initialize_hwe  (bool probs_fixed);

    /// Initializes my_parameters for the nhwe option
    ///
    void  initialize_nhwe (bool probs_fixed);
    //@}
    
    // Data members.  
    sm_option my_option;
    double    my_freq_A;                 // Frequency of allele A.
    double    my_probs[NUM_OF_PROBS];    // Genotype probabilities.
    bool      my_default;                // True if user doesn't specify a value for any estimable parameter.
};


} 
} 
#include "segreg/freq_sub_model.ipp"

#endif


