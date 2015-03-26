#ifndef SEGREG_TYPE_SUB_MODEL_H
#define SEGREG_TYPE_SUB_MODEL_H
//============================================================================
// File:      type_sub_model.h
//                                                                          
// Author:    Geoff Wedig (wedig@darwin.cwru.edu)
//                                                                          
// History:   0.1 gcw Initial Implementation                             May 2001
//                Baechle   reformatted, continued development           Jun 2001
//                djb renamed, factored out genotype_sub_model  
//                    and added genotype_specific_variance_sub_model.    7-30-01
//                djb added genotype_specific_susceptibility_sub_model.  1-16-02
//                                                                          
// Notes:     defines the genotype specific sub-models for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <math.h>
#include "error/internal_error.h"
#include "segreg/sub_model_base.h"
#include "segreg/definitions.h"

namespace SAGE
{

namespace SEGREG {

extern const std::string  TYPE_MEAN_NAME;
extern const std::string  TYPE_SUSCEPT_NAME;
extern const std::string  TYPE_VAR_NAME;

/// \name genotype specific mean sub-model constants
///
/// There are no boundaries on either means or suscepts per Yi Dong
/// 1-17-02 (ie, boundaries are +/- infinity)

//@{
extern const double MEAN_DEFAULT_VALUE;
extern const bool   MEAN_DEFAULT_FIXED;
//@}

/// \name genotype specific variance sub-model constants
//@{
extern const double VAR_DEFAULT_VALUE;
extern const bool   VAR_DEFAULT_FIXED;
extern const double VAR_EPSILON;
extern const double VAR_LB;
//@}

// We don't want default constructors, so:
//lint -esym(1712, genotype_specific_*sub_model)

//----------------------------------------------------------------------------
//  Class:    genotype_specific_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
/// base class for all genotype specific sub_models.
/// This class should never be directly instantiated.

class genotype_specific_sub_model : public SegregSubmodel
{
  public:

    friend class genotype_specific_mean_susc_sub_model;
    friend class genotype_specific_mean_sub_model;
    friend class genotype_specific_susceptibility_sub_model;
    friend class genotype_specific_variance_sub_model;
    friend class CovariateSubmodel;

    /// Options for the derived sub_models.  
    /// Note that some options are only available for certain sub_models. 
    /// Specifically, three_dec and three_inc are only valid for mean and
    /// susceptibility models.

    enum sm_option 
    {
      one = 1,       ///< One mean/var/susc.
      two,           ///< Two means/vars/suscs
      three,         ///< Three means/vars/suscs
      two_dom,       ///< Two means/vars/suscs with uAB = uAA
      two_rec,       ///< Two means/vars/suscs with uAB = uBB
      three_add,     ///< Three means/vars/suscs with ...
      three_dec,     ///< Three means/vars/suscs with uAA > uAB > uBB
      three_inc      ///< Three means/vars/suscs with uAA < uAB < uBB
    };

    static string  option_2_parameter  (sm_option option);

    // Gets.
    sm_option option() const;
    string    option_description() const;  
    string    name() const;

//    double  get_parameter(size_t gt) const;
    double  parameter(genotype_index gt) const;

    /// Returns the current number of types (1, 2 or 3)
    ///
    size_t get_type_count() const;
    
    /// Control the management of the two option.  By default, two acts like
    /// a dominant model, but can be set otherwise if needed.
    //@{

    bool  is_two_dom() const;
    //@}

    /// Returns \c true if the current option is any of two, two_dom,
    /// or two_rec; \c false otherwise.
    bool is_two_option() const;

    /// Returns \c true if the current option is any of three, three_add,
    /// three_dec, or three_inc; \c false otherwise.
    bool is_three_option() const;
    
    bool is_default() const { return my_default; }

    /// Tests the submodel for completeness, which is defined as all
    /// parameters specified as finite values.
    bool is_complete() const;

  protected:

    enum derived_type { mean, variance, susc, invalid };

    // Constructors protected to prevent construction of base class

    genotype_specific_sub_model(derived_type type_word, cerrorstream& errors = sage_cerr);
    genotype_specific_sub_model(const genotype_specific_sub_model& other);
    genotype_specific_sub_model&  operator=(const genotype_specific_sub_model& other);

    static string  static_option_description(sm_option option, derived_type type);

    /// \brief Set options
    ///
    /// These options set the model to some predefined status, returning
    /// true if successful and false otherwise.
    //@{
    bool  set_one(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);
    bool  set_two(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);
    bool  set_two_dom(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);
    bool  set_two_rec(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);

    bool  set_three(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);
    bool  set_three_add(const model_input& type_AA, const model_input& type_AB, const model_input& type_BB);
    bool  set_three_dec(const model_input& mean_AA, const model_input& mean_AB, const model_input& mean_BB);
    bool  set_three_inc(const model_input& mean_AA, const model_input& mean_AB, const model_input& mean_BB);

    void  set_defaults();
    bool  set_one_one_value(const model_input& mi);
    bool  set_one_two_values(const model_input& mi_one, const model_input& mi_two);
    void  set_three_add_set_AB();

    /// sets three values for three, three_dec and three_inc models.
    /// pstatus is the status used if all are not fixed.
    bool  set_three(const model_input& type_AA, const model_input& type_AB,
                    const model_input& type_BB, bool non_fixed_status);
    //@}

    /// \brief Maximization Initialization
    ///
    //@{
      
    /// This function synchronizes the my_parameters vector with what's in
    /// my_types and my_types fixed.  Necessary for inclusion in maxfun
    /// or other times where that must be synched.
    virtual int finalizeConfiguration();

    /// Initialize my_parameters for one type
    ///
    void  initialize_one();

    /// Initialize my_parameters for two types with A dominant over B (\f$t_AB = t_AA\f$) 
    ///
    void  initialize_two_dom();

    /// Initialize my_parameters for two types with A recessive under B (\f$t_AB = t_BB\f$) 
    ///
    void  initialize_two_rec();

    /// Initialize my_parameters for three types 
    ///
    void  initialize_three();

    /// Initialize my_parameters for three_add option (\f$t_AB = (t_AA + t_BB) / 2\f$)
    ///
    void  initialize_three_add();
    //@}

    /// \brief Option specific synchronization functions (sync w. Maxfun).
    ///
    /// Note that not all of these are used for every model.
    //@{
    int  synchronize_one       ();
    int  synchronize_two       ();
    int  synchronize_three     ();
    int  synchronize_two_dom   ();
    int  synchronize_two_rec   ();
    int  synchronize_three_add ();
    int  synchronize_three_dec ();
    int  synchronize_three_inc ();
    //@}    

    /// Various type to string options, short or long variants in singular and plural.
    //@{
    string  short_sing_type  () const;
    string  long_sing_type   () const;
    string  short_plural_type() const;
    string  long_plural_type () const;
    string  type_param       () const;

    static string  short_sing_type  (derived_type t);
    static string  long_sing_type   (derived_type t);
    static string  short_plural_type(derived_type t);
    static string  long_plural_type (derived_type t);
    static string  type_param       (derived_type t);
    //@}

    /// Constraint Checking
    //@{
    typedef bool (genotype_specific_sub_model::* constraint_met_ptr)() const;

    bool  constraint_met(constraint_met_ptr constraint);
    bool  three_dec_constraint_met() const;
    bool  three_inc_constraint_met() const;
    //@}
        
    // Data members.  
    mutable derived_type  my_derived_type; ///< Enum defining if it's mean, var or susc

    sm_option  my_option;
    int        my_type_count;

    double     my_types[NUM_OF_TYPES];
    bool       my_types_fixed[NUM_OF_TYPES];

    bool  my_default;                  ///< True if user gives no value for any estimable parameter.
    bool  my_two_is_dom;               ///< Option two is interpreted as two dominant, otherwise interpretted as two_rec
};

//----------------------------------------------------------------------------
//  Class:    genotype_specific_mean_susc_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
/// base class for mean and susceptibility models
class genotype_specific_mean_susc_sub_model : public genotype_specific_sub_model
{
  public:

    genotype_specific_mean_susc_sub_model&  operator=(const genotype_specific_mean_susc_sub_model& other);

    virtual ~genotype_specific_mean_susc_sub_model();
    
    void  dump(std::ostream& out) const;
    
    // Sets.
    bool  set(sm_option opt, const model_input& mean_AA, const model_input& mean_AB, 
              const model_input& mean_BB);
    bool  set(sm_option opt, const model_input& mean_AA, const model_input& mean_AB, 
              const model_input& mean_BB, primary_type pt);
              
  protected:

    // Constructor/destructor.  Protected to avoid accidental creation.
    genotype_specific_mean_susc_sub_model(derived_type t, cerrorstream& errors = sage_cerr);
    genotype_specific_mean_susc_sub_model(const genotype_specific_mean_susc_sub_model& other);

    virtual int update();
  
    // Ancillary functions.
    void  type_fixity(bool types[]) const;
    bool  primary_trait_type_ok(primary_type pt) const;

  private:
  
    friend class genotype_specific_variance_sub_model;

//    vector<ParameterInput>& get_params();
};

//----------------------------------------------------------------------------
//  Class:    genotype_specific_mean_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class genotype_specific_mean_sub_model : public genotype_specific_mean_susc_sub_model
{
  public:
  
    // Constructor/destructor.  
    genotype_specific_mean_sub_model(cerrorstream& errors = sage_cerr);
    genotype_specific_mean_sub_model(const genotype_specific_mean_sub_model& other);
    genotype_specific_mean_sub_model&  operator=(const genotype_specific_mean_sub_model& other);
    virtual ~genotype_specific_mean_sub_model();
    
    static string  option_2_description(sm_option option);
};


//----------------------------------------------------------------------------
//  Class:    genotype_specific_susceptibility_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class genotype_specific_susceptibility_sub_model : public genotype_specific_mean_susc_sub_model
{
  public:
                     
    // Constructor/destructor.  
    genotype_specific_susceptibility_sub_model(cerrorstream& errors = sage_cerr);
    genotype_specific_susceptibility_sub_model(const genotype_specific_susceptibility_sub_model& other);
    genotype_specific_susceptibility_sub_model&  operator=(const genotype_specific_susceptibility_sub_model& other);
    virtual ~genotype_specific_susceptibility_sub_model();
    
    static std::string  option_2_description(sm_option option);
};


//----------------------------------------------------------------------------
//  Class:    genotype_specific_variance_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class genotype_specific_variance_sub_model : public genotype_specific_sub_model
{
  public:
    friend class model;
  
    // Constructor/destructor.  
    genotype_specific_variance_sub_model(genotype_specific_mean_susc_sub_model* m_ptr, 
                                         cerrorstream& errors = sage_cerr);
    genotype_specific_variance_sub_model(const genotype_specific_variance_sub_model& other);
    genotype_specific_variance_sub_model&  operator=(const genotype_specific_variance_sub_model& other);
    virtual ~genotype_specific_variance_sub_model();
    
    // Gets.
    
    void  dump(std::ostream& out) const;
    

    // Sets.
    bool  set(sm_option opt, const model_input& var_AA, const model_input& var_AB, 
              const model_input& var_BB, bool type_missing, bool trans_missing);
    bool  set(sm_option opt, const model_input& var_AA, const model_input& var_AB, 
              const model_input& var_BB, bool type_missing, bool trans_missing, primary_type pt);
              
    static std::string  option_2_description(sm_option option);

    // Mean/variance consistency.
    bool  mv_types(sm_option var_option) const;
    
  protected:
    
    virtual int update();
    
    // Ancillary functions.
    void  mv_types_message(sm_option opt, sm_option mean_opt);
    bool  primary_trait_type_ok(primary_type pt) const;
    bool  var_meets_constraints(model_input& var);
  
    // Data members.  
    genotype_specific_mean_susc_sub_model*  my_m_ptr;
};


} 
} 

#include "segreg/type_sub_model.ipp"

#endif
