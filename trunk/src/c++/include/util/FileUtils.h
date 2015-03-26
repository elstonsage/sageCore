#ifndef UTIL_FILE_UTILS_H
#define UTIL_FILE_UTILS_H

#include <iostream>
#include <string>
#include <fstream>

namespace SAGE {
namespace UTIL {

/// Cool file utility functions.
class FileUtils
{
public:

  ///
  /// Returns \c true if the given file exists, \c false otherwise.
  /// \param filename The name of the file
  static bool fileExists(const std::string & filename);

  ///
  /// Reads from the ifstream until it hits a newline or EOF.
  /// \param i The ifstream from which to read
  /// \returns The string read in (trailing newline is stripped of course!)
  static std::string getLine(std::ifstream & i);
};

inline bool
FileUtils::fileExists(const std::string & filename)
{
  std::ifstream f;
  
  f.open(filename.c_str(), std::ios::in);
  
  return f.is_open();
}

inline std::string
FileUtils::getLine(std::ifstream & i)
{
  std::string s = "";
  char        c;
  
  while(i.get(c))  
  {
    if(c != '\n')
      s += c;
    else
      break;
  }
                            
  return s;
}

} // End namespace UTIL
} // End namespace SAGE

#endif
