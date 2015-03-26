#ifndef DATATYPES_H
#define DATATYPES_H

#include "error/internal_error.h"
#include "output/Output.h"
#include "maxfun/maxfun.h"

namespace SAGE   {
/// MAXFUN namespace
namespace MAXFUN {

// Forward class declarations:

class SequenceCfg;
class Datatypes;
class DebugCfg;
class ParameterMgr;
class APIMaxFunction;
class Maximizer;
class Parameter;
class OutputFormatter;
class LogitTransformer;
class Transformer;
class ParameterIterator;
class ParameterConstIterator;
class ParamCalculator;
class Results;
class Submodel;
class NewSubmodel;

const string GLOBAL_GROUP_NAME = "GLOBAL";

/// \hideinitializer MF_INFINITY is available for use as a parameter's upper or lower bound.
const double MF_INFINITY = std::numeric_limits<double>::infinity();

/// @name Conversion functions
//@{

  ///
  /// Converts a parameter name (composed of a group name and a parameter name) into a
  /// single string. The two names are separated by a colon.
  /// \param group_name The group to which the parameter belongs.
  /// \param param_name The name of the parameter.
  string paramname2str(const string & group_name, const string & param_name);

  ///
  /// Converts a string in which group and individual names are colon-delimited into a
  /// pair of strings (the first string is group name, the second is individual parameter name).
  /// \param combined_name The name that will be split into group and individual names.
  pair<string, string> str2paramname(const string & combined_name);

//@}

} // End namespace MAXFUN
} // End namespace SAGE

#endif
