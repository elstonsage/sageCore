#ifndef LODLINK_MLE_SUB_MODEL_H
#define LODLINK_MLE_SUB_MODEL_H
//============================================================================
// File:      mle_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   09/04/02 - created.                                   djb
//            10/25/06 - modified to use a different Maxfun 
//                       interface.                                 djb
//                                                                          
// Notes:     defines the mle sub-model for LODLINK.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "maxfunapi/maxfunapi.h"

namespace SAGE
{

namespace LODLINK
{

typedef MAXFUN::Parameter::ParamTypeEnum  pstatus;

extern const string MLE_NAME;

extern const string  AVERAGE_THETA_DESC;
extern const pstatus  THETA_DEFAULT_STATUS;
extern const string  MALE_THETA_DESC;
extern const string  FEMALE_THETA_DESC;
extern const double  THETA_INIT_VALUE;
extern const double  THETA_LOWER_BOUND;
extern const double  THETA_STRICT_UPPER_BOUND;
extern const double  THETA_UPPER_BOUND;

extern const string  ALPHA_DESC;
extern const pstatus  ALPHA_DEFAULT_STATUS;
extern const double  ALPHA_INIT_VALUE;
extern const double  ALPHA_LOWER_BOUND;
extern const double  ALPHA_UPPER_BOUND;

enum one_theta_idx { AVERAGE, ALPHA_ONE };
enum two_theta_idx { MALE, FEMALE, ALPHA_TWO };
enum sex { male, female };

//----------------------------------------------------------------------------
//  Class:    mle_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class mle_sub_model : public MAXFUN::Submodel
{
  public:
    
    friend std::ostream& operator<<(std::ostream& out, const mle_sub_model& mle);   
    
    // Constructor/destructor.  
    mle_sub_model(bool ss = false, bool ua = false, cerrorstream& errors = sage_cerr);
    
    // Gets.
    virtual string  option_description() const;  
    virtual string  name() const;
    const std::vector<MAXFUN::ParameterInput>&  parameters() const;
    double          average_theta() const;
    double          average_theta_ub() const;
    double          male_theta() const;
    double          male_theta_ub() const;
    double          female_theta() const;
    double          female_theta_ub() const;
    double          theta(LODLINK::sex s) const;
    double          alpha() const;
    bool            is_sex_specific() const;
    bool            uses_alpha() const;
    const MAXFUN::ParameterMgr*  parameter_mgr() const;
    
    // Sets.
    void  set(bool ss = false, bool ua = false);
    void  reset();
    void  set_strict_limits();
    void  set_relaxed_limits();
    void  fix_alpha();
    void  unfix_alpha();
    void  constrain_thetas();
    void  unconstrain_thetas();
    void  set_average_theta(double value);
    void  set_male_theta(double value);
    void  set_female_theta(double value);
    void  set_alpha(double value);
    
  protected:
    int update();
    
  private:
    void init();
    void  internally_synchronize();
    
    // Data members.  
    double  my_average_theta;     // theta is recombination fraction.
    double  my_male_theta;
    double  my_female_theta;
    double  my_alpha;             // alpha is proportion of families w. linkage.
    
    bool  sex_specific;
    bool  use_alpha;
};

#include "mle_sub_model.ipp"
} 
} 

#endif

