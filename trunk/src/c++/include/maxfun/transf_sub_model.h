#ifndef TRANSF_SUB_MODEL_H
#define TRANSF_SUB_MODEL_H
//============================================================================
// File:      transformation_sub_model.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/2/01 - created.                                   djb
//                                                                          
// Notes:     defines the transformation sub-model for SEGREG.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "error/internal_error.h"
#include "LSF/LSF.h"
#include "app/aparser.h"
#include "sub_model.h"

namespace SAGE
{

/// @name transformation sub-model constants
//@{
extern const std::string  TRANSFORMATION_NAME;
extern const double       LAMBDA_ONE_DEFAULT_VALUE;
extern const bool         LAMBDA_ONE_DEFAULT_FIXED;
extern const double       LAMBDA_ONE_DEFAULT_LB;
extern const double       LAMBDA_ONE_DEFAULT_UB;
extern const double       LAMBDA_TWO_DEFAULT_VALUE;
extern const bool         LAMBDA_TWO_DEFAULT_FIXED;
//@}

//----------------------------------------------------------------------------
//  Class:    transformation_sub_model
//                                                                          
//----------------------------------------------------------------------------
//
class transformation_sub_model : public sub_model
{
  public:
    enum sm_option { box_cox = 1, george_elston, no_trans };
    
    // Constructor/destructor.  
    inline transformation_sub_model(cerrorstream& errors = sage_cerr);
    inline transformation_sub_model(const transformation_sub_model& other);
    inline transformation_sub_model&  operator=(const transformation_sub_model& other);
    virtual inline ~transformation_sub_model();
    
    // Gets.
    inline sm_option       option() const;
    virtual inline string  option_description() const;  
    virtual inline string  name() const;
    inline double          lambda_one() const;
    inline double          lambda_two() const;
    
    inline void  dump(std::ostream& out) const;

    // Calculation.
    inline bool  calculate_geom_mean(const std::vector<double>& traits) const;
    
    inline bool  transform(std::vector<double>& traits) const;
    inline bool  transform(double&              trait ) const;

    // Sets.
    bool  set(sm_option opt, const model_input& lambda_one, const model_input& lambda_two, 
              double lambda_one_lb, double lambda_one_ub);

    inline void  set_option(sm_option opt);

    inline void  set_clear_each_sync(bool clear=true) const { clear_each_sync = clear; }
                  
    static inline std::string  option_2_description(sm_option option);
    static inline std::string  option_2_parameter(sm_option option);
    
  protected:
    virtual int synchronize(parameter_iterator start);
    inline void  internally_synchronize();
    
    // Ancillary functions.
    inline bool  set_none(const model_input& lambda_one, const model_input& lambda_two,
                   double lambda_one_lb, double lambda_one_ub);
    inline bool  set_lambda_one_limits(const model_input& lambda_one, double lambda_one_lb, double lambda_one_ub);
    inline bool  set_lambda_one(const model_input& lambda_one);
    inline bool  set_lambda_two(const model_input& lambda_two);
    inline void  adjust_lambda_one();
                   
    // Transformation.
    static inline double  sign(double value);
    double  bc_geom_mean(const std::vector<double>& traits) const;
    double  ge_geom_mean(const std::vector<double>& traits) const;
    inline void  bc_transform_power_zero     (double G, double& traits) const;
    inline void  bc_transform_power_non_zero (double G, double& traits) const;
    inline void  ge_transform_power_zero     (double G, double& traits) const;
    inline void  ge_transform_power_non_zero (double G, double& traits) const;
    
    // Data members.  
    sm_option  my_option;
    double     my_lambda_one;
    double     my_lambda_two;

    mutable bool clear_each_sync;

    // For calculation purposes
    
    mutable double my_geometric_mean;
};

bool  parse_transformation_sub_model
  (transformation_sub_model&           tsm,
   const LSFBase*                      param,
   transformation_sub_model::sm_option option = transformation_sub_model::box_cox,
   cerrorstream&                       errors = sage_cerr);

} 

#include "maxfun/transf_sub_model.ipp"

#endif
