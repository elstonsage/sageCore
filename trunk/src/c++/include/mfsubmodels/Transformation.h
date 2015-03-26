// Added pvalue means to configuration class.   -djb 2/23/9


#ifndef MFSUBMODELS_TRANSFORMATION_H
#define MFSUBMODELS_TRANSFORMATION_H

#include <cmath> 
#include "output/Output.h"
#include "numerics/log_double.h"
#include "maxfunapi/maxfunapi.h"
#include "error/internal_error.h"
#include "app/aparser.h"
//#include "app/output_streams.h"

namespace SAGE        {
namespace MFSUBMODELS {

/// \brief The transformation toolkit wrapper class
///
/// The transformation toolkit provides transformation of traits using
/// one of three models: 'none' (no transformation), box-cox or george-elston.
///
/// For the latter two models, there are two parameters, lambda1 and lambda2,
/// which are used in the transformation.  By default, lambda1 is
/// unfixed and varies from -1.0 to infinity (initial estimate 1.0), while 
/// lambda2 is fixed and varies from negative infinity to infinity (initial estimate 0.0).
class Transformation
{
  public:
  
  /// @name Enums
  //@{
    
    /// \brief Describes which transformation option is chosen.
    enum TransformationType
    {
      BOX_COX = 1,          ///< Box-Cox Transformation
      GEORGE_ELSTON,        ///< George-Elston Transformation
      NONE                  ///< No transformation
    }; 
      
  //@}
  
  class Parser;

  /// \brief Describes the configuration options for a transformation block.
  ///
  /// This configuration object contains information on the model selected as
  /// well as the lambda1 and lambda2 parameters.  Access to the lambda1 and 2
  /// parameters is limited to changin the fixed/unfixed status and the 
  /// initial values of the parameters.
  class Configuration
  {
    friend class Parser;
    
    public:
    
    /// @name Constructors / operators
    //@{
    
      Configuration();
      Configuration(const Configuration & other);

      Configuration& operator= (const Configuration & other);
      
    //@}

    /// @name Getting
    //@{
      TransformationType get_type() const;
      const MAXFUN::ParameterInput& get_lambda1() const;
      const MAXFUN::ParameterInput& get_lambda2() const;
      OUTPUT::Table summarize_as_table() const;
      OUTPUT::Section dump() const;

    //@}   
      
    /// @name Setting
    //@{
    
      void set_type(TransformationType t);
      
      bool set_lambda1_init_est(double d);
      bool set_lambda1_fixed(bool is_fixed = false);
      void set_lambda1_lower_bound(double d);
      void set_lambda1_upper_bound(double d);

      bool set_lambda2_init_est(double d);
      bool set_lambda2_fixed(bool is_fixed = false);
      void  set_lambda1_pvalue_mean(double pvalue_mean);
      void  set_lambda2_pvalue_mean(double pvalue_mean);
      double  get_lambda1_pvalue_mean() const;
      double  get_lambda2_pvalue_mean() const;
      bool  either_lambda_estimated() const;
      bool  both_lambdas_estimated() const;
            
      
    //@}   
    private:
      TransformationType     my_type;
      MAXFUN::ParameterInput my_lambda1;
      MAXFUN::ParameterInput my_lambda2;
      double  my_lambda1_pvalue_mean;
      double  my_lambda2_pvalue_mean;
  };
  
  /// Parses an LSFBase* into a Configuration object.
  class Parser
  {
    public:
      static Configuration parse_block(const LSFBase* param, std::ostream & info, cerrorstream & errors);
        static void parse_option(Configuration& config, const LSFBase* param, std::ostream & info, cerrorstream & errors);
        static void parse_lambda1(Configuration& config, const LSFBase* param, std::ostream & info, cerrorstream & errors);
        static void parse_lambda2(Configuration& config, const LSFBase* param, std::ostream & info, cerrorstream & errors);
  };
  
  /// \brief Provides an interface to the transformation functionality.
  ///
  /// The transformation's facade provides access to the transformation functions.
  /// There are two phases to using these functions.  First, the geometric
  /// mean must be calculated.  This must be done before values are transformed.
  /// After that, values or groups of values can be transformed.
  /// (Q: Clear cached geometric mean option added to the maxfun function, maybe?)
  class Facade
  {
    public:
    
      Facade(const Configuration& config, MAXFUN::Function & func);
      Facade(const Facade& other);
    
      bool calculate_geometric_mean(const vector<double>& traits, bool calculate = true);
      bool transform(vector<double>& traits) const;
      bool transform(double& trait) const;
      const Configuration& get_configuration() const;
      double  geometric_mean() const;

    private:
    
      // Disallowed
      Facade& operator= (const Facade& other);

      typedef bool (Facade::*CalcFcn)(const vector<double>& traits);
      typedef bool (Facade::*TransfFcn)(double& trait)  const;
    
      bool calculate_geom_mean_none(const vector<double>& traits);
      bool calculate_geom_mean_box_cox(const vector<double>& traits);
      bool calculate_geom_mean_george_elston(const vector<double>& traits);
    
      bool transform_none(double& trait) const;
      bool transform_box_cox(double& trait) const;
      bool transform_george_elston(double& trait) const;
    
      double get_lambda1() const;
      double get_lambda2() const;
    
      TransformationType     my_type;
      CalcFcn                my_calc_fcn;
      TransfFcn              my_transf_fcn;
      MAXFUN::ParameterMgr&  my_mgr;
      size_t                 my_lambda1_id;
      size_t                 my_lambda2_id;
      
      double                 my_geom_mean;
      Configuration          my_config;
  };
};

} // End namespace MFSUBMODELS
} // End namespace SAGE

#include "mfsubmodels/Transformation.ipp"

#endif

