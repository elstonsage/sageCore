#ifndef SEGREG_ASCER_SUB_MODEL_H
#define SEGREG_ASCER_SUB_MODEL_H
//============================================================================
// File:      ascertainment_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/18/01  - created.                          djb
//                                                                          
// Notes:     defines the ascertainment sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "sub_model_base.h"
#include "fped/fped.h"
#include "app/output_streams.h"

namespace SAGE
{

namespace SEGREG
{

/// @name ascertainment sub-model constants
//@{
extern const std::string  ASCERTAINMENT_NAME;
extern const double       ASCER_INCLUDE_DEFAULT_VALUE;
extern const double       ASCER_INDIC_THRESH_DEFAULT_VALUE;
//@}

//----------------------------------------------------------------------------
//  Class:    ascertainment_sub_model
//                                               
//  Note: **** onset option not implemented
//                           
//----------------------------------------------------------------------------
//
class ascertainment_sub_model
{
  public:
    enum s_sm_option  { none, founders, psf, founders_plus_psf };
    
    // - 'not_specified' is for internal use, not an actual option.
    //
    enum v_sm_option  { not_specified, actual, gte_thresh, lte_thresh, thresh_indic,
                        onset };
    
    // Constructor/destructor.  
    ascertainment_sub_model(cerrorstream& errors = sage_cerr);
    ascertainment_sub_model(const ascertainment_sub_model& other);
    ascertainment_sub_model&  operator=(const ascertainment_sub_model& other);
    virtual ~ascertainment_sub_model();
    
    // Gets.
    s_sm_option                 s_option() const;
    v_sm_option                 v_option() const;

    virtual string              option_description()   const;
    string                      s_option_description() const;
    string                      v_option_description() const;

    virtual string              name() const;
    string                      psf_indicator() const;
    const std::vector<double>&  includes() const;
    string                      thresh_indicator() const;
    double                      thresh() const;
    double                      thresh_high() const;
    double                      thresh_low() const;
    double                      indic_thresh() const;

    void  dump(std::ostream& out) const;

    /// Tests to determine if an individual is included in our ascertained set
    //@{
    /// determines if an individual is in the C set.  Note that this test
    /// does not test for individual validity with regards to the primary
    /// trait and/or covariates. This tests based upon the s_sm_option type. 
    /// When type is none, returns false.
    bool is_ind_in_C (FPED::MemberConstPointer m) const;

    /// Tests an individual for being in the proband sampling frame
    bool is_ind_in_PSF       (FPED::MemberConstPointer m) const;

    /** Returns the type of the individual.
     *  This type is one of four states:
     *    - not_specified (missing),
     *    - actual,
     *    - gte_thresh,
     *    - lte_thresh.
     *  This fuction first determines if the individual is in C (if not,
     *  type is missing) then determines if it's in PSF (if not, it's a
     *  founder, and therefore uses actual) then determines which threshold
     *  type, if any to use. */
    v_sm_option get_ind_type      (FPED::MemberConstPointer m) const;

    //@}
    
    /** Calculates the high and low thresholds, if needed.  If there's a problem,
     *  (can't calculate, high < low) then it returns false. */
    bool calculate_thresholds(const FPED::Multipedigree& rmp, size_t primary_trait);
    
    // Sets.
    bool  set(s_sm_option so, v_sm_option vo,
              const RPED::MultiPedigree* mp, const string& pi, const std::vector<double>& includes, 
              const string& ti, double thresh, double thresh_high, double thresh_low, double it); 
    

    static string  s_option_2_description(s_sm_option option);
    static string  s_option_2_parameter(s_sm_option option);
    static string  v_option_2_description(v_sm_option option);
    static string  v_option_2_parameter(v_sm_option option);
  
  private:

    bool calculate_threshold_high(const FPED::Multipedigree& rmp, size_t primary_trait);
    bool calculate_threshold_low (const FPED::Multipedigree& rmp, size_t primary_trait);

    /** If thresh_indic is used, determines if the member has a valid entry.
     *  Returns true if thresh_indic is not used or it is used and 
     *  the individual has a valid entry. */
    bool   does_ind_have_valid_thresh     (FPED::MemberConstPointer m) const;

    double get_ind_thresh                 (FPED::MemberConstPointer m) const;
    double get_ind_trait_value            (FPED::MemberConstPointer m,
                                           const string& trait) const;

    bool ind_uses_thresh_high(FPED::MemberConstPointer m) const;
    bool ind_uses_thresh_low (FPED::MemberConstPointer m) const;  

    // Set functions.
    bool  set_none(v_sm_option vo, const RPED::MultiPedigree* mp, const string& pi, const std::vector<double>&, 
                   const string& ti, double thresh, double thresh_high, double thresh_low, double it);
    bool  set_founders(v_sm_option vo, const RPED::MultiPedigree* mp, const string& pi, const std::vector<double>&, 
                       const string& ti, double thresh, double thresh_high, double thresh_low, double it);
    bool  set_psf(v_sm_option vo, const RPED::MultiPedigree* mp, const string& pi, const std::vector<double>&, 
                  const string& ti, double thresh, double thresh_high, double thresh_low, double it);
    bool  set_founders_plus_psf(v_sm_option vo, const RPED::MultiPedigree* mp, const string& pi, const std::vector<double>&, 
          const string& ti, double thresh, double thresh_high, double thresh_low, double it);
    
    // Ancillary functions.
    bool  set_actual(const RPED::MultiPedigree* mp, const string& ti, 
                     double thresh, double thresh_high, double thresh_low, double low);
    bool  set_gte_thresh(const RPED::MultiPedigree* mp, const string& ti, 
                         double thresh, double thresh_high, double thresh_low, double it);
    bool  set_lte_thresh(const RPED::MultiPedigree* mp, const string& ti, 
                         double thresh, double thresh_high, double thresh_low, double it);
    bool  set_thresh_indic(const RPED::MultiPedigree* mp, const string& ti, 
                           double thresh, double thresh_high, double thresh_low, double it);
    void  set_includes(const std::vector<double>& includes);
    bool  trait_valid(const string& trait_name, RPED::RefTraitInfo::trait_t req_type1, 
                      RPED::RefTraitInfo::trait_t req_type2, const RPED::MultiPedigree* mp);
    
    // Data members. 
    
    /// Error stream to which we send errors about configuration options.
    ///
    cerrorstream my_errors;
    
    s_sm_option  my_s_option;          // Condition set option.
    v_sm_option  my_v_option;          // Condition value option.
    string  my_psf_indicator;          // Proband sampling frame indicator.
    std::vector<double>  my_includes;          // Used for set_option, psf_indic.
    string  my_thresh_indicator;       // Threshold indicator (continuous trait).
    mutable double  my_thresh;                 // Used for val_options, thresh_high and thresh_low.
    mutable double  my_thresh_high;            // Used for val_option, thresh_indic.
    mutable double  my_thresh_low;             // Used for val_option, thresh_indic.
    mutable double  my_indic_thresh;           // Used for val_option, thresh_indic.
    bool  my_t_specified;              // Did user supply value for thresh?
    bool  my_th_specified;             // Did user supply value for thresh_high?
    bool  my_tl_specified;             // Did user supply value for thresh_low?
};

inline ostream& operator<< (ostream& out, const ascertainment_sub_model& sm);

}
}

#include "segreg/ascertainment_sub_model.ipp"

#endif


