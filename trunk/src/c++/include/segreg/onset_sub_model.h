#ifndef SEGREG_ONSET_SUB_MODEL_H
#define SEGREG_ONSET_SUB_MODEL_H
//============================================================================
// File:      onset_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/28/01  - created.                          djb
//                                                                          
// Notes:     defines the age of onset sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "error/internal_error.h"
#include "sub_model_base.h"

namespace SAGE
{

namespace SEGREG
{

extern const std::string  ONSET_NAME;

//----------------------------------------------------------------------------
//  Class:    onset_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class onset_sub_model
{
  public:
    enum type_option 
    {
      t_A,  ///< Age   is major type dependent
      t_S   ///< Susc. is major type dependent
    };
    
    enum multi_option
    {
      m_N, ///< Neither is polygenically dependent 
      m_A, ///< Age     is polygenically dependent 
      m_S  ///< Susc.   is polygenically dependent 
    };

    enum ageon_option
    {
      a_AN, ///< Age   dep on MT, Neither dep on poly.
      a_AA, ///< Age   dep on MT and poly.
      a_AS, ///< Age   dep on MT, Susc.   dep on poly.
      a_SN, ///< Susc. dep on MT, Neither dep on poly.
      a_SA, ///< Susc. dep on MT, Age     dep on poly.
      a_SS  ///< Susc. dep on MT and poly. 
    };
    
    // Constructor/destructor.  
    onset_sub_model(cerrorstream& errors = sage_cerr);
    onset_sub_model(const onset_sub_model& other);
    onset_sub_model&  operator=(const onset_sub_model& other);
    virtual ~onset_sub_model();
    
    // Gets.
    type_option     t_option() const;
    multi_option    m_option() const;
    ageon_option    a_option() const;
    string          t_option_description() const;
    string          m_option_description() const;
    string          option_description() const;
    string          name() const;
    string          affection_status() const;
    string          age_of_onset() const;
    string          age_at_exam() const;
    
    void  dump(std::ostream& out) const;
    
    // Sets.
    bool  set(type_option t_option, multi_option m_option, const RPED::RefMultiPedigree* mp,  
              const string& status, const string& onset, const string& exam); 
  
    static string  t_option_2_description(type_option option);
    static string  t_option_2_parameter(type_option option);
    static string  m_option_2_description(multi_option option);
    static string  m_option_2_parameter(multi_option option);
    bool           no_poly_loci(); // due to JA
    static  bool   static_no_poly_loci(); // due to JA
  
  protected:
    
    // Ancillary functions.
    bool  trait_valid(const string& trait_name, 
                      RPED::RefTraitInfo::trait_t req_type, const RPED::RefMultiPedigree* mp);
    bool  trait_given(const string& trait_name, const string& parameter);

    /// Sets the a option according to the t and m options
    void  set_a_option();

    // Data members. 
    mutable cerrorstream my_errors;

    type_option   my_t_option;
    multi_option  my_m_option;
    ageon_option  my_a_option;
    string     my_affection_status;    // Name of a binary trait.
    string     my_age_of_onset;        // Name of a continuous trait.
    string     my_age_at_exam;         // Name of a continuous trait.
};

inline ostream& operator<< (ostream& out, const onset_sub_model& sm);

} 
} 

#include "segreg/onset_sub_model.ipp"

#endif


