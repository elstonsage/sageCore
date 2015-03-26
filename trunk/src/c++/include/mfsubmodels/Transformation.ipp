#ifndef MFSUBMODELS_TRANSFORMATION_H
#include "mfsubmodels/Transformation.h"
#endif

namespace SAGE        {
namespace MFSUBMODELS {

inline
Transformation::Configuration::Configuration()
{
  my_type    = NONE;
  my_lambda1 = MAXFUN::ParameterInput("Transformation", "Lambda1", MAXFUN::Parameter::INDEPENDENT, 1.0, -1.0,                 MAXFUN::MF_INFINITY);
  my_lambda2 = MAXFUN::ParameterInput("Transformation", "Lambda2", MAXFUN::Parameter::FIXED,       0.0, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);
  my_lambda1_pvalue_mean = 0;
  my_lambda2_pvalue_mean = 0;
}

inline
Transformation::Configuration::Configuration(const Configuration & other) :
  my_type    (other.my_type),
  my_lambda1 (other.my_lambda1),
  my_lambda2 (other.my_lambda2),
  my_lambda1_pvalue_mean(other.my_lambda1_pvalue_mean),
  my_lambda2_pvalue_mean(other.my_lambda2_pvalue_mean)
{}

inline
Transformation::Configuration& Transformation::Configuration::operator=
    (const Configuration & other)
{
  if(this != &other)
  {
    my_type    = other.my_type;
    my_lambda1 = other.my_lambda1;
    my_lambda2 = other.my_lambda2;
    my_lambda1_pvalue_mean = other.my_lambda1_pvalue_mean;
    my_lambda2_pvalue_mean = other.my_lambda2_pvalue_mean;
  }
  
  return *this;
}

inline
Transformation::TransformationType
    Transformation::Configuration::get_type() const 
{
  return my_type;
}

inline
const MAXFUN::ParameterInput&
  Transformation::Configuration::get_lambda1() const
{
  return my_lambda1;
}

inline
const MAXFUN::ParameterInput&
  Transformation::Configuration::get_lambda2() const
{
  return my_lambda2;
}

inline
void Transformation::Configuration::set_type(TransformationType t)
{
  my_type = t;
}

inline
bool Transformation::Configuration::set_lambda1_init_est
  (double d)
{
  my_lambda1.initial_estimate = d;

  return true;
}

inline bool  
Transformation::Configuration::either_lambda_estimated() const
{
  return  my_type != Transformation::NONE     &&
          (my_lambda1.initial_type == MAXFUN::Parameter::INDEPENDENT ||
           my_lambda2.initial_type == MAXFUN::Parameter::INDEPENDENT   );
}

inline bool  
Transformation::Configuration::both_lambdas_estimated() const
{
  return  my_type != Transformation::NONE                           &&
          my_lambda1.initial_type == MAXFUN::Parameter::INDEPENDENT &&
          my_lambda2.initial_type == MAXFUN::Parameter::INDEPENDENT   ;
}

inline
bool Transformation::Configuration::set_lambda1_fixed
  (bool is_fixed)
{
  if(is_fixed) my_lambda1.initial_type = MAXFUN::Parameter::FIXED;
  else         my_lambda1.initial_type = MAXFUN::Parameter::INDEPENDENT;

  return true;
}

inline
void Transformation::Configuration::set_lambda1_lower_bound(double d)
{
  my_lambda1.lower_bound = d;
}

inline
void Transformation::Configuration::set_lambda1_upper_bound(double d)
{
  my_lambda1.upper_bound = d;
}            

inline
bool Transformation::Configuration::set_lambda2_init_est (double d)
{
  my_lambda2.initial_estimate = d;

  return true;
}

inline
bool Transformation::Configuration::set_lambda2_fixed(bool is_fixed)
{
  if(is_fixed) my_lambda2.initial_type = MAXFUN::Parameter::FIXED;
  else         my_lambda2.initial_type = MAXFUN::Parameter::INDEPENDENT;

  return true;
}

inline void  
Transformation::Configuration::set_lambda1_pvalue_mean(double pvalue_mean)
{
  my_lambda1_pvalue_mean = pvalue_mean;
}

inline void  
Transformation::Configuration::set_lambda2_pvalue_mean(double pvalue_mean)
{
  my_lambda2_pvalue_mean = pvalue_mean;
}

inline double  
Transformation::Configuration::get_lambda1_pvalue_mean() const
{
  return  my_lambda1_pvalue_mean;
}

inline double  
Transformation::Configuration::get_lambda2_pvalue_mean() const
{
  return  my_lambda2_pvalue_mean;
}            

inline
Transformation::Facade::Facade
    (const Configuration& config, MAXFUN::Function& func)
  : my_type(config.get_type()), 
    my_mgr(func.getMgr()),
    my_lambda1_id(static_cast<size_t>(-1)),
    my_lambda2_id(static_cast<size_t>(-1)),
    my_geom_mean(0.0),
    my_config(config)
{
  if(my_type != NONE)
  {
    MAXFUN::ParameterInput lambda1 = config.get_lambda1();
    MAXFUN::ParameterInput lambda2 = config.get_lambda2();
    
    my_lambda1_id = my_mgr.addParameter(lambda1);
    my_lambda2_id = my_mgr.addParameter(lambda2);
    
    my_mgr.getParameter(my_lambda1_id).setPValueInf(MAXFUN::Parameter::TWO_SIDED, config.get_lambda1_pvalue_mean());
    my_mgr.getParameter(my_lambda2_id).setPValueInf(MAXFUN::Parameter::TWO_SIDED, config.get_lambda2_pvalue_mean());
  }
  
  switch(my_type)
  {
    case NONE :          my_calc_fcn   = &Facade::calculate_geom_mean_none;
                         my_transf_fcn = &Facade::transform_none;
                         break;
    case BOX_COX :       my_calc_fcn   = &Facade::calculate_geom_mean_box_cox;
                         my_transf_fcn = &Facade::transform_box_cox;
                         break;
    case GEORGE_ELSTON : my_calc_fcn   = &Facade::calculate_geom_mean_george_elston;
                         my_transf_fcn = &Facade::transform_george_elston;
                         break;
  }
}

inline
Transformation::Facade::Facade(const Facade & other)
  : my_type       (other.my_type),
    my_calc_fcn   (other.my_calc_fcn),
    my_transf_fcn (other.my_transf_fcn),
    my_mgr        (other.my_mgr),
    my_lambda1_id (other.my_lambda1_id),
    my_lambda2_id (other.my_lambda2_id),
    my_geom_mean  (other.my_geom_mean),
    my_config(other.my_config) 
{ }

inline
const Transformation::Configuration&
Transformation::Facade::get_configuration() const
{
  return  my_config;
}

inline
double
Transformation::Facade::get_lambda1() const
{
  return my_mgr(my_lambda1_id);
}

inline
double
Transformation::Facade::get_lambda2() const
{
  return my_mgr(my_lambda2_id);
}
    
}
}

