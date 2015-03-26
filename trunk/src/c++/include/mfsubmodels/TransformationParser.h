#ifndef MFSUBMODELS_TRANSFORMATIONPARSER_H
#define MFSUBMODELS_TRANSFORMATIONPARSER_H

#include "mfsubmodels/TransformationSubmodel.h"

namespace SAGE        {
namespace MFSUBMODELS {

/** \brief Parses a TransformationSubmodel
  *
  */
class TransformationParser
{
  public:
 
    ///
    /// Parses a TransformationSubmodel.
    /// \param tsm The TransformationSubmodel instance to set up
    /// \param param The LSFBase* pointing to the transformation sub-block
    /// \param option The transformation method (optional, defaults to box_cox)
    /// \param errors The errorstream to use (optional, defaults to sage_cerr)
    static bool parse(
      TransformationSubmodel& tsm,
      const LSFBase* param,
      TransformationSubmodel::sm_option option = TransformationSubmodel::box_cox,
      cerrorstream&                     errors = sage_cerr);
      
    ///
    /// Parses a NewTransformationSubmodel.
    /// \param tsm The TransformationSubmodel instance to set up
    /// \param param The LSFBase* pointing to the transformation sub-block
    /// \param option The transformation method (optional, defaults to box_cox)
    /// \param errors The errorstream to use (optional, defaults to sage_cerr)
    static bool parse(
      NewTransformationSubmodel& tsm,
      const LSFBase* param,
      NewTransformationSubmodel::sm_option option = NewTransformationSubmodel::box_cox,
      cerrorstream&                     errors = sage_cerr);
  private:
  
    static bool get_value(
      model_input&   mi,
      double         def,
      const LSFBase* param,
      const string&  name_phrase,
      cerrorstream&  errors);

    static void  get_lambda_one_bounds(
      double&        lambda_one_lb,
      double&        lambda_one_ub,
      const LSFBase* param, 
      const string&  name_phrase,
      bool           fixed,
      cerrorstream&  errors);
                 
    static bool  get_lambda_one(
      model_input&   lambda_one,
      double&        lambda_one_lb,
      double&        lambda_one_ub,
      const LSFBase* param,  
      cerrorstream&  errors);
                               
    static bool  get_lambda_two(
      model_input&   lambda_two,
      const LSFBase* param,
      cerrorstream&  errors);
};

} // End namespace MFSUBMODELS
} // End namespace SAGE

#endif
