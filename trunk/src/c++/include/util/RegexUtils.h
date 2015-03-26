#ifndef UTIL_REGEX_UTILS_H
#define UTIL_REGEX_UTILS_H

#include "boost/regex.hpp"
#include <vector>
#include <string>
#include <iostream>

namespace SAGE {
namespace UTIL {

class RegexUtils
{
public:

  ///
  /// Searches the given input string for the match pattern.
  /// \param pattern The regular expression pattern to search for
  /// \param input_line The string to search
  /// \returns \c true if the pattern is found, \c false otherwise
  static bool doesPatternMatch(
    const std::string & pattern,
    const std::string & input_line);

  ///
  /// Searches the given input string for the match pattern.
  /// \param pattern The regular expression pattern to search for
  /// \param input_line The string to search
  /// \returns A vector of found patterns; the first element is the entire found pattern;
  /// each subsequent element is the found subexpression. If the pattern string is malformed
  /// or no results are found, returns an empty vector.
  static std::vector<std::string> matchPattern(
    const std::string & pattern, 
    const std::string & input_line);

  ///
  /// Searches the given input string for the match pattern.
  /// This function is basically identical to matchPattern(), except
  /// that it will only return the whole matched string. It is useful
  /// if you don't have any subexpressions in your pattern.
  /// \param pattern The regular expression pattern to search for
  /// \param input_line The string to search
  /// \returns The whole matched string, or an empty string if the pattern
  /// was invalid or the pattern didn't match.
  static std::string matchPatternOnce(
    const std::string & pattern, 
    const std::string & input_line);

  static std::string searchAndReplace(
    const std::string & input_line,
    const std::string & search_pattern,
    const std::string & replacement_text);
};

inline bool 
RegexUtils::doesPatternMatch(
    const std::string & pattern,
    const std::string & input_line)
{
  return matchPattern(pattern, input_line).size();
}

inline std::vector<std::string>
RegexUtils::matchPattern(
    const std::string & pattern,
    const std::string & input_line)
{
  // Set up variables:
  
  std::vector<std::string> results (0);
  boost::cmatch            matches;
  boost::regex             regex_pattern(pattern);

  // Match pattern & stick results in the results vector:
  
  try
  {
    if(boost::regex_match(input_line.c_str(), matches, regex_pattern) == true)
      for(boost::cmatch::const_iterator i = matches.begin(); i != matches.end(); ++i)
        results.push_back(*i);
  }
  catch(const std::exception & e) 
  { 
    results.clear();
  }

  // Return results:

  return results;
}

inline std::string 
RegexUtils::matchPatternOnce(
    const std::string & pattern, 
    const std::string & input_line)
{
  std::vector<std::string> results = matchPattern(pattern, input_line);
  
  if(results.size() == 0)
    return "";
  else
    return results[0];
}

inline std::string 
RegexUtils::searchAndReplace(
    const std::string & input_line,
    const std::string & search_pattern,
    const std::string & replacement_text)
{
  return regex_replace(input_line, boost::regex(search_pattern), replacement_text);
}


} // End namespace UTIL
} // End namespace SAGE

#endif
