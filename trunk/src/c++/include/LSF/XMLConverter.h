#ifndef LSF_XML_CONVERTER_H
#define LSF_XML_CONVERTER_H

#include <string>
#include <sstream>
#include "LSF/LSF.h"

namespace SAGE {
namespace LSF  {

/// \brief Converts an LSFBase * into an XML-compliant string
///
///  
class XMLConverter
{
public:
  ///
  /// Given an LSFBase pointer, returns its entire contents as an XML-formatted string.
  static std::string toXML(const LSFBase * data);

private:

  ///
  /// Prepends a backslash to characters that SHOULD have backslashes accordinging to the XML
  /// standard.
  static std::string convertAttributeVal(const std::string & txt);
};

} // End namespace LSF
} // End namespace SAGE

#endif
