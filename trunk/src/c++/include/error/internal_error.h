#ifndef INTERNAL_ERROR_H
#define INTERNAL_ERROR_H

#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"

namespace SAGE {

/// @name Internal error functions
//@{

/// Produces a fatal warning message then exits the program.
/// \param o The outputstream to which output will be directed.
/// \param file Name of the file in which the error was generated.
/// \param line The line on which the error was generated.
inline void internal_error(cerrorstream& o, const std::string& file, int line, std::string txt)
{
  o << priority(fatal) << "Sorry, S.A.G.E. has detected an unexpected error detected at ("
    << file << '/' << line << ") and cannot continue.  "
    << "Especially if maximization is involved, possible reasons for this include: "
    << "\n1) A large number of the trait values are identical.  "
    << "Examine the distribution of your data."
    << "\n2) A specified transformation cannot be performed.  "
    << "Check initial values, if specified, and/or for outliers and try again."
    << "\n3) S.A.G.E. is not trying appropriate initial parameter estimates.  "
    << "If you have reason to believe you have good initial estimates for some or all the parameters, try them."
    << std::endl
    << txt << std::endl;

  exit(1);
}

/// Produces a fatal warning message then exits the program.
/// \param file Name of the file in which the error was generated.
/// \param line The line on which the error was generated.
inline void internal_error(const std::string& file, int line, std::string txt)
{
  internal_error(sage_cerr, file, line, txt);
}

#define SAGE_internal_error() SAGE::internal_error(__FILE__, __LINE__, "")
#define SAGE_internal_error_msg(txt) SAGE::internal_error(__FILE__, __LINE__, txt)

#define SAGE_DBG std::cout << "DEBUG (" << __FILE__ << ") Reached line #" << __LINE__ << std::endl;
#define SAGE_PAUSE SAGE_DBG; { char t; std::cout << "<pause>"; std::cin >> t; } SAGE_DBG;

//@}

} // End namespace SAGE

#endif
