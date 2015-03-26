#ifndef PARAMETERINPUT_H
#define PARAMETERINPUT_H
//===========================================================
//
//  File:  ParameterInput.h
//
//  Author:  Stephen Gross
//
//  Copyright 2004 R. C. Elston
//===========================================================


#include "output/Output.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/Parameter.h"

namespace SAGE   {
namespace MAXFUN {

/** \class ParameterInput
 *  \brief A helper class for adding new parameters.
 *
 */
class ParameterInput
{
  public:
  /// @name Constructors & operators
  //@{
    /// Default constructor
    ParameterInput();

    /// Alternate constructor
    ParameterInput(string    group_name, 
                   string    param_name, 
                   Parameter::ParamTypeEnum initial_type, 
                   double    initial_estimate, 
                   double    lower_bound, 
                   double    upper_bound);

    /// Copy constructor
    ParameterInput(const ParameterInput &);

    /// Assignment operator
    ParameterInput& operator=(const ParameterInput &);

    // Copy operation
    void copy(const ParameterInput &);

  //@}

  /// @name Basic information
  //@{

    /// The name of the group to which the parameter belongs.
    std::string group_name;

    /// The name of the parameter.
    std::string param_name;

    /// The initial estimate of the parameter (default is 0).
    double initial_estimate;

    /// The initial type of the parameter (default is MAXFUN::Parameter::INDEPENDENT).
    Parameter::ParamTypeEnum initial_type;

    /// The lower bound of the parameter.
    double lower_bound;

    /// The upper bound of the parameter.
    double upper_bound;

    /// Lookup index for parameter value.
    /// When the maximization is actively taking place, you can use the index value
    /// as a lookup index for the parameter's value:
    /// \code getMaxfunInfo()->operator()(index) \endcode will return a double reference
    /// to the parameter's current estimate.
    int index;
    
  //@}
  
  /// @name Debugging
  //@{
  
    OUTPUT::Table dump() const
    {
      return (OUTPUT::Table("ParameterInput")
        << (OUTPUT::TableRow() << "group_name"       << group_name)
        << (OUTPUT::TableRow() << "param_name"       << param_name)
        << (OUTPUT::TableRow() << "initial_estimate" << initial_estimate)
        << (OUTPUT::TableRow() << "initial_type"     << ParamTypeEnum2str(initial_type))
        << (OUTPUT::TableRow() << "lower_bound"      << lower_bound)
        << (OUTPUT::TableRow() << "upper_bound"      << upper_bound)
        << (OUTPUT::TableRow() << "index"            << index));
    }
  
  //@}
};

} // End namespace MAXFUN
} // End namespace SAGE

#endif
