#ifndef SEGREG_PARSER_H
#define SEGREG_PARSER_H
//============================================================================
// File:      parser.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/14/01 created                     - djb
//                                                                          
// Notes      Defines class, parser, for parsing segreg_analysis block of
//            a SAGE parameter file.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include "rped/rped.h"
#include "app/aparser.h"
#include "LSF/LSF.h"
#include "LSF/Attr.h"
#include "segreg/sub_model_base.h"
#include "segreg/type_sub_model.h"
#include "segreg/freq_sub_model.h"
#include "segreg/resid_sub_model.h"
#include "segreg/transmission_sub_model.h"
#include "segreg/cov_sub_model.h"
#include "segreg/model.h"
#include "segreg/fpmm_sub_model.h"
#include "segreg/onset_sub_model.h"
#include "segreg/ascertainment_sub_model.h"

namespace SAGE
{

namespace SEGREG
{

//----------------------------------------------------------------------------
//  Class:    parser
//                                                                          
//  Purpose:  extract information from segreg_analysis block of SAGE parameter
//            file.
//                                                                          
//----------------------------------------------------------------------------
//
class parser : public SAGE::APP::BasicParser
{
  public:
    
    struct base_ptrs
    {
      base_ptrs();
    
      const LSFBase*  title;
      const LSFBase*  trait;
      const LSFBase*  composite_trait;
      const LSFBase*  type_mean;
      const LSFBase*  type_var;
      const LSFBase*  type_suscept;
      const LSFBase*  mean_cov;
      const LSFBase*  var_cov;
      const LSFBase*  suscept_cov;
      const LSFBase*  class_;              // 'class' a reserved word.
      const LSFBase*  fpmm;
      const LSFBase*  resid;
      const LSFBase*  transformation;
      const LSFBase*  geno_freq;
      const LSFBase*  transmission;
      const LSFBase*  ascertainment;
      const LSFBase*  prev_constraint;
      const LSFBase*  prev_constraints;   // For multiple constraits
      const LSFBase*  prev_estimate;
      const LSFBase*  output_options;
      const LSFBase*  output;
      
      const LSFBase*  maxfun_options;     ///< Stores maxfun details.
    };
    
    typedef CovariateSubmodel::CovariateTypeEnum covariate_type;

  
    // Constructor/destructor.
    parser(const RPED::RefMultiPedigree* mp, ostream& messages = std::cout,
           SAGE::cerrorstream& errors = SAGE::sage_cerr);
  
    virtual ~parser() { my_mp = 0; }

    virtual void  parse_symbols(const SymbolTable* syms);
    virtual void  parse_parameter(const LSFBase* param);
    virtual void  parse_test_parameter_section(const LSFBase* params);
    virtual void  parse_test_parameter(const LSFBase* param);
    
    const model& get_model() const;  
    long   analysis_number() const;

    bool        zero_poly_loci; // due to JA (needed for no polygenic loci in fpmm)
    bool        type_var_one; // due to JA (needed in polygenic loci in fpmm)
    double      like_cutoff; // due to JA, for tuning acceptance of competing starting values
  private:
    void  init_parse();
    void  print_title();
    void  print_footer();
    void  parameter_duplicated(const std::string& param_name);

    /// classify parameter by type (one of the base_ptrs), and set the base
    /// ptr.  If already set, model is invalid.
    void  classify_parameter(const LSFBase* param);

    void  set_parameter(const string& name, const LSFBase** ptr, const LSFBase* param);

    void  second_pass_parse();

    // General ancillary functions.
    void  reset();
    void  init_genotype_map();
    bool  get_attributes(model_input& mi, double def, const LSFBase* param, 
                         const string& name_phrase, bool value_open) const;
    RPED::RefTraitInfo::trait_t  primary_trait_type() const;
    void  not_relevant_error_message(const string& sub_model_name, primary_type type);
                         
    // Type (mean, var and suscept) sub-block ancillary functions.
    bool  get_types(model_input types[], const LSFBase* param, 
                    const std::string& keyword);                  // Uses non-const member func., parse_string().     
                                                                 
    // Freq sub-block ancillary functions.
    void  get_probs_fixed(bool& probs_fixed, bool& probs_fixed_specified, const LSFBase* param);
    void  get_probs(double probs[], const LSFBase* param);
    void  get_freq_A(double& freq_A, const LSFBase* param);
    
    // Covariate sub-block ancillary functions.
    bool  get_covariate(model_input& coeff, bool& interaction, 
                                covariate_type type, const LSFBase* param); //  const; changed in making interaction
                                                                           // not available
    void  reset_covariate_sub_model(covariate_type type);
    void  get_interaction(bool& interaction, const LSFBase* param, 
                                covariate_type type, const string& name_phrase); //  const;

    // Transmission sub-block ancillary functions.
    bool  get_taus(model_input taus[], const LSFBase* param);
    void  get_no_bounds(bool& no_bounds, const LSFBase* param);
    
    // Prevalence sub-block ancillary functions.
    double  get_real_attribute(const LSFBase* param, const std::string& attribute,
                                                     const std::string& parameter, 
                                                     const std::string& sub_block);
    size_t  get_int_attribute(const LSFBase* param, const std::string& attribute,
                                                    const std::string& parameter, 
                                                    const std::string& sub_block);
    bool  prev_attribute_invalid(double value, const std::string& name);
    bool  R_N_constraint_not_met(double R, double N);
    
    // Model class ancillary functions.
    void  type_warning(const string& type);
        
    // Parsers.
    void  parse_mean_sub_model(const LSFBase* param, const std::string& keyword);
    void  parse_var_sub_model(const LSFBase* param);
    void  parse_freq_sub_model(const LSFBase* param);
    void  parse_resid_sub_model(const LSFBase* param);
    void  parse_transmission_sub_model(const LSFBase* param);
    void  parse_transformation_sub_model(const LSFBase* param);
    void  parse_cov_sub_model(const LSFBase* param, covariate_type type);
    void  parse_fpmm_sub_model(const LSFBase* param);
    void  parse_onset_sub_model(const LSFBase* param);
    void  parse_ascertainment_sub_model(const LSFBase* param);
    void  parse_prev_constraint(const LSFBase* param);
    void  parse_prev_constraints(const LSFBase* param);
    void  parse_prev_estimate_sub_model(const LSFBase* param);    
    void  parse_type_suscept_sub_model(const LSFBase* param);

    void  parse_primary_trait(const LSFBase* param);
    void  parse_title(const LSFBase* param);
    void  parse_output_parameter(const LSFBase* param);
    void  parse_output_attribute(const LSFBase* param);
    void  parse_model_class(const LSFBase* param);
    void  parse_output_options(const LSFBase* param);

    // Some simple helper functions
    
    inline bool model_valid() const;
    inline bool process_block( const LSFBase* ptr) const;
    
    // Meta-constraints(sub-block/sublock or extra sub-block/sublock requirements).
    void  check_meta_constraints();
    void  no_primary_trait();
    void  print_covariate_non_exclusive_error(const CovariateSubmodel& sm1,
                                              const CovariateSubmodel& sm2) const;

    void  check_ignore_3();
    void  check_ignore_4();
    void  check_ignore_6();
    
    // Data members.
    long idnum; // Analysis count

    base_ptrs ptrs;

    const RPED::RefMultiPedigree*  my_mp;
    model  my_model;

    bool  transmission_parsed;                   // Used to insure geno_freq not parsed after transmission.
    bool  fpmm_parsed;
    bool  resid_parsed;
    bool  transformation_parsed;
    
    bool  file_name_root_parsed;

    std::map<string, unsigned>  my_genotype_map;

    // Output
    
    ostream& messages;

    //lint --e{1712}
};

std::ostream&  operator<<(std::ostream& out, const parser::base_ptrs& ptrs);
  

}
}

#include "segreg/parser.ipp"

#endif

