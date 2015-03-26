#ifndef APP_VERSION_NUMBER_H
#define APP_VERSION_NUMBER_H

#include <string>
#include <sstream>

namespace SAGE {
namespace APP  {

/// \brief Represents a version of SAGE.
///
/// \par Overview
///
/// A version is composed of four numbers: main version, sub version,
/// micro version, and beta version (where a beta version of 0 indicates
/// that the version is \b not a beta version).
///
/// Examples include: 1.2.3.4, 5.0.1.0, 2.2.3.119, and so on.
///
/// \par Install path
///
/// There is a standard way to convert a version number into an install path.
/// This is because SAGE is installed in a standard location according to its
/// version. See VersionNumber::toInstallPath() for more details.
class VersionNumber
{
public:

  ///
  /// Constructor.
  VersionNumber(int main, int sub, int micro, int beta);  

  ///
  /// Renders this version as an install path
  /// (eg: (non-beta) "SAGE 1.2.3", (beta) "SAGE 3.2.4 Beta 4")
  std::string toInstallPath() const;

  ///
  /// Renders this version in pretty print.
  /// (eg: (non-beta) "S.A.G.E. v1.2.3", (beta) "S.A.G.E. 3.2.4 Beta 4")
  std::string toPrettyPrint() const;

  /// Indicates whether or not this is a beta version (beta version != 0, that is).
  bool isBeta() const;

  /// Returns the main version.
  int getMainVersion() const;
  
  /// Returns the sub version.
  int getSubVersion() const;
 
  /// Returns the micro version.
  int getMicroVersion() const;
  
  /// Returns the beta version (0 = not beta)
  int getBetaVersion() const;

private:

  int my_main_version;
  int my_sub_version;
  int my_micro_version;
  int my_beta_version;
};

extern VersionNumber version;

inline int VersionNumber::getMainVersion  () const { return my_main_version;  }
inline int VersionNumber::getSubVersion   () const { return my_sub_version;   }
inline int VersionNumber::getMicroVersion () const { return my_micro_version; }
inline int VersionNumber::getBetaVersion  () const { return my_beta_version;  }

inline
VersionNumber::VersionNumber(int main, int sub, int micro, int beta)
{
  my_main_version  = main;
  my_sub_version   = sub;
  my_micro_version = micro;
  my_beta_version  = beta;
}

inline bool VersionNumber::isBeta() const { return my_beta_version != 0; }

inline std::string
VersionNumber::toInstallPath() const
{
  std::ostringstream s;
  
  s << "SAGE " 
    << my_main_version << "." 
    << my_sub_version << "."
    << my_micro_version;
    
  if(isBeta())
    s << " Beta " << my_beta_version;
    
  return s.str();
}

inline std::string
VersionNumber::toPrettyPrint() const
{
  std::ostringstream s;
  
  s << "S.A.G.E. v" 
    << my_main_version << "." 
    << my_sub_version << "."
    << my_micro_version;
    
  if(isBeta())
    s << " Beta " << my_beta_version;
    
  return s.str();
}

} // End namespace APP
} // End namespace SAGE


#endif
