#ifndef TRANSFORMATIONSUBMODEL_H
#define TRANSFORMATIONSUBMODEL_H
//============================================================================
// File:      TransformationSubmodel.h
//                                                                          
// Author:    Dan Baechle
//            Stephen Gross
//                                                                          
// History:   7/2/01 - created.                                   djb
//            5/25/04 - adapted to new MaxfunAPI                  sag
//                                                                          
// Notes:     defines the transformation sub-model for SEGREG.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "maxfunapi/Submodel.h"
#include "error/internal_error.h"
#include "app/aparser.h"

namespace SAGE   {
namespace MAXFUN {

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

/** \class TransformationSubmodel
 *  \brief Transforms stuff.
 *
 *  Given a vector of doubles, transform each member of the vector
 *  per D21 of Yi Dong's write-up dated 8/3/01.
 *
 *  For Box-Cox option trait value + lambda two must be greater
 *  than 0 for every trait or return value is false and no trans-
 *  formations are performed.
 *
 *  Note that although the model for the log transform is specified such that
 *  lambda1 is equal to 0, for numerical reasons, we use the log form of the
 *  equations any time we're within a small range around 0.  This avoids
 *  numerical instabilities that occur when multiplying and dividing numbers that
 *  are very close to 0.
 */
class TransformationSubmodel : public Submodel
{
public:

  enum sm_option { box_cox = 1, george_elston, no_trans };
    
  /// @name Constructors & operators
  //@{

    ///
    /// Constructor
    inline TransformationSubmodel(cerrorstream& errors = sage_cerr);

    ///
    /// Copy constructor
    inline TransformationSubmodel(const TransformationSubmodel & other);

    ///
    /// Assignment operator
    inline TransformationSubmodel&  operator=(const TransformationSubmodel & other);

    //
    // Copy operation (for copy constructor and operator=)
    inline void copy(const TransformationSubmodel & other);

    ///
    /// Destructor
    virtual inline ~TransformationSubmodel();

  //@}
    
  /// \internal
  /// @name Virtual interface
  //@{

    /// \internal
    /// Updates any necessary values for a new function evaluation.
    virtual inline int update();
    
  //@}

  /// @name Accessors
  //@{
    ///
    /// Returns the type of transformation being used.
    inline sm_option option() const;

    ///
    /// Returns a string describing the type of transformation being used.
    virtual inline string option_description() const;  

    ///
    /// Returns the name of the submodel (in this case, "TRANSFORMATION")
    virtual inline string name() const;

    ///
    /// Returns the current value of lambda 1.
    inline double lambda_one() const;

    ///
    /// Returns the current value of lambda 2.
    inline double lambda_two() const;
    
    ///
    /// Dumps runtime debugging info to the indicated ostream.
    /// \param out The ostream to which output should be directed.
    inline void dump(std::ostream& out) const;

  //@}

  /// @name Calculation functions
  //@{

    ///
    /// Calculates a geometric mean.
    /// \param traits The vector of values on the basis of which the geometric mean will be calculated.
    /// \retval true Mean was successfully calculated.
    /// \retval false Mean was \b not successfully calculated.
    inline bool calculate_geom_mean(const std::vector<double> & traits) const;
    
    ///
    /// Transforms a vector of values.
    /// \param traits The vector of values which will be transformed.
    /// \retval true Transformation was successful.
    /// \retval false Transformation was \b not successful.
    bool transform(std::vector<double>& traits) const;

    ///
    /// Transforms a single value.
    /// \param trait The value which will be transformed.
    /// \retval true Transformation was successful.
    /// \retval false Transformation was \b not successful.
    bool transform(double & trait) const;

  //@}

  /// @name Mutators
  //@{

    ///
    /// Sets up the basic transformation information.
    /// \param opt Which transformation type to use.
    /// \param lambda_one Default value for lambda1.
    /// \param lambda_two Default value for lambda2.
    /// \param lambda_one_lb Lower bound for lambda1.
    /// \param lambda_one_ub Upper bound for lambda1.
    bool set(  sm_option opt, 
        const model_input& lambda_one, 
        const model_input& lambda_two, 
        double lambda_one_lb, 
        double lambda_one_ub);

    ///
    /// Sets the basic transformation option.
    /// \param opt Which transformation to use.
    void set_option(sm_option opt);

    void set_clear_each_sync(bool clear=true) const { clear_each_sync = clear; }

  //@}
                  
  /// @name Static helper functions
  //@{

    static std::string  option_2_description(sm_option option);

    static std::string  option_2_parameter(sm_option option);

      //@}

protected:
  // Ancillary functions.

    bool  set_none(const model_input& lambda_one, const model_input& lambda_two,
                   double lambda_one_lb, double lambda_one_ub);
    bool  set_lambda_one_limits(const model_input& lambda_one, double lambda_one_lb, double lambda_one_ub);
    bool  set_lambda_one(const model_input& lambda_one);
    bool  set_lambda_two(const model_input& lambda_two);
    void  adjust_lambda_one();
    
    bool lambda_one_close_to_zero() const;
                   
  // Transformation

    static double  sign(double value);
    double  bc_geom_mean(const std::vector<double>& traits) const;
    double  ge_geom_mean(const std::vector<double>& traits) const;
    void  bc_transform_power_zero     (double G, double& traits) const;
    void  bc_transform_power_non_zero (double G, double& traits) const;
    void  ge_transform_power_zero     (double G, double& traits) const;
    void  ge_transform_power_non_zero (double G, double& traits) const;

  virtual int finalizeConfiguration();
    
  // Data members

    sm_option         my_option;
    mutable bool      clear_each_sync;
    mutable double    my_geometric_mean;
    
    double my_lambda_one;
    double my_lambda_two;
};

bool parseTransformationSubmodel(
         TransformationSubmodel&           tsm,
   const LSFBase*                          param,
         TransformationSubmodel::sm_option option = TransformationSubmodel::box_cox,
         cerrorstream&                     errors = sage_cerr);

} // End namespace MAXFUN
} // End namespace SAGE

#include "maxfunapi/TransformationSubmodel.ipp"

#endif
