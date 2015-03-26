#ifndef UTIL_STRING_UTILS_H
#define UTIL_STRING_UTILS_H

//===================================================================
//
//  Standard includes
//
//===================================================================

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>

#include "error/internal_error.h"
#include "util/RegexUtils.h"

namespace SAGE {
namespace UTIL {

//-----------------------------------------------------------------------
//  Original author:    Paul J. Weiss (http://www.codeproject.com/string/stringsplit.asp)
//------------------------------------------------------------------------

/// \brief Provides functions for manipulating strings (including date conversion)
class StringUtils
{
  public:

  /// @name Case conversion
  //@{
  
    ///
    /// Converts the given string to uppercase.
    static std::string toUpper(const std::string & input);
  
  //@}

  /// @name String tokenizing
  //@{

    ///
    /// Splits the given string on the given delimiter. Places the results in the
    /// given string vector. If no instances of the delimiter are found, places
    /// the entire input string in the results vector.
    /// \param input The source string
    /// \param delimiter The string that delimits the occurences you're searching for.
    /// \param results The string vector in which the found results will be placed
    /// \returns An int describing the number of instances found
    static size_t splitString(const std::string& input, const std::string& delimiter, std::vector<std::string>& results);

    ///
    /// Splits the given string on the given list of single-character delimiters. Places the results in the
    /// given string vector. If no instances of the delimiter are found, places
    /// the entire input string in the results vector.
    ///
    /// Please note that a straight sequence of delimiter characters encountered in the input string will be treated as
    /// a \b single delimiter. For instance, consider the following code:
    ///
    /// \code
    /// std::string input = "this,is,a;,;,string";
    /// std::vector<std::string> results;
    /// size_t i = splitMultiDelimitedString(input, ",;", results);
    /// \endcode
    ///
    /// The integer 'i' will equal four, because the ';,;,' sequence between 'a' and 'string' will be interpreted as a 
    /// single delimiting instance.
    /// \param input The source string
    /// \param delimiter The string in which each character is considered a possible delimiter; if set to an empty string (""),
    /// the function will split the string on whitespace.
    /// \param results The string vector in which the found results will be placed
    /// \returns An int describing the number of instances found
    static size_t splitMultiDelimitedString(const std::string& input, const std::string& delimiter, std::vector<std::string>& results);

    ///
    /// Splits the given string on all whitespace. This is equivalent to splitMultiDelimitedString with 
    /// \c delimiter = " \t\n".
    /// \param input The source string
    /// \param results The string vector in which the found results will be placed
    /// \returns An int describing the number of instances found
    static size_t splitOnWhitespace(const std::string & input, std::vector<std::string> & results);

    ///
    /// Indents every line of the given input string by 'depth' number of spaces.
    /// \param src The input string
    /// \param depth The number of spaces to indent each line.
    /// \param skip_first_line If \c true, this function will NOT indent the first line of the text; if false (the default), it will.
    static std::string indentString(const std::string & src, size_t depth, bool skip_first_line = false);
    
    ///
    /// Reformats a string so that (1) all contiguous whitespace is replaced with single-space, and (2) no more than
    /// \c width characters occur between newlines. If width < 1, this function will treat width as 1.
    static std::string lineWrapString(const std::string & src, size_t width);

  //@}

  /// @name String/time conversion
  //@{

    ///
    /// Converts a time_t to a string.
    /// \param t The time_t to convert
    /// \returns The string form of the given time_t
    static std::string convertTimeToStr(const time_t & t);

    ///
    /// Converts a string to a time_t.
    /// The string must take the form: "Thu Dec 16 15:05:41 2004".
    ///
    /// Note #1: The month is read in case-insensitive, so you can have any
    /// kind of capitalization in the month that you want.
    ///
    /// Note #2: If you don't want to specify the day (Mon-Fri), the day string
    /// can be expressed as '---'. For example: "--- Dec 16 15:05:41 2004" 
    ///   
    /// \param t The string to convert
    /// \param time_inst The time_t to place the converted values into
    /// \retval 0 Conversion successful
    /// \retval 1 Conversion not successful (usually due to a malformed string)
    static int convertStrToTime(const std::string & t, time_t & time_inst);
    
  //@}
  
  /// @name Numeric conversion
  //@{
  
    ///
    /// If a string represents a number in scientific format, its exponent may have been
    /// rendered with three digits (ie: 1.25e-009). If you want to strip the leading zero
    /// from that exponent, use this function:
    static std::string stripLeadingExpZero(const std::string & s);
  
  //@}
  
};

inline std::string 
StringUtils::toUpper(const std::string & input)
{
  std::string s = input;
  
  for(size_t i = 0; i < input.length(); ++i)
    s[i] = toupper(input[i]);
    
  return s;
}

//===========================================
//
//  splitMultiDelimitedString(...)
//
//===========================================
inline size_t 
StringUtils::splitMultiDelimitedString(const std::string& input, const std::string& delimiter_list, std::vector<std::string>& results)
{
  // If the delimiter is empty, just split on whitespace:
  
  if(delimiter_list == "")
    return splitString(input, delimiter_list, results);

  // Replace all delimiters with the last delimiter:

  std::string final_delimiter = delimiter_list.substr(delimiter_list.length() - 1, 1);
  std::string adj_input       = "";

  // Replace multi-delimiter sequences with single-delimiter sequences of the final delimiter:

  for(size_t i = 0; i < input.length(); ++i)
  {
    char cur_char = input[i];

    // Have we hit a delimiter?
    if(delimiter_list.find(cur_char) != std::string::npos)
    {
      adj_input += final_delimiter;
      
      // Advance past this (possibly multi-character) delimiting sequence:
      
      while(delimiter_list.find(input[i]) != std::string::npos && i < input.length())
        i++;
        
      i--;
    }
    else
    {
      adj_input += cur_char;
    }
  }

  return splitString(adj_input, final_delimiter, results);
}

//===========================================
//
//  splitString(...)
//
//===========================================
inline size_t 
StringUtils::splitString(const std::string& input, const std::string& delimiter, std::vector<std::string>& results)
{
  results.clear();

  if(delimiter == "")
    return 0;
  
  size_t cur_idx = 0;
  
  while(cur_idx < input.length())
  {
    size_t next_delimiter_idx = input.find(delimiter, cur_idx),
           add_till_idx       = next_delimiter_idx == std::string::npos ? input.length() : next_delimiter_idx;

    if(add_till_idx - cur_idx)
      results.push_back(input.substr(cur_idx, add_till_idx - cur_idx));

    cur_idx = add_till_idx + delimiter.length();
  }

  return results.size();
}

//=======================================
//
//  splitOnWhitespace(...)
//
//=======================================
inline size_t 
StringUtils::splitOnWhitespace(const std::string & input, std::vector<std::string> & results)
{
  return splitMultiDelimitedString(input, " \t\n", results);
}


//=============================================================
//
//  StringUtils::convertTimeToStr(...)
//
//=============================================================
inline std::string 
StringUtils::convertTimeToStr(const time_t & t)
{
  std::string t_str = ctime(&t);

  if(t_str[t_str.length() - 1] == '\n')
    t_str = t_str.substr(0, t_str.length() - 1);

  return t_str;
}

//=============================================================
//
//  StringUtils::convertStrToTime(...)
//
//=============================================================
inline int
StringUtils::convertStrToTime(const std::string & t, time_t & new_time)
{
//            10        20
//  012345678901234567890123456
// [Thu Dec 16 15:05:41 2004]

  if(t.length() != 24)
    return 1;

  tm temp_tm;

  // Process month:

  std::string month = StringUtils::toUpper(t.substr(4,3));

       if(month == "JAN") temp_tm.tm_mon = 0;
  else if(month == "FEB") temp_tm.tm_mon = 1;
  else if(month == "MAR") temp_tm.tm_mon = 2;
  else if(month == "APR") temp_tm.tm_mon = 3;
  else if(month == "MAY") temp_tm.tm_mon = 4;
  else if(month == "JUN") temp_tm.tm_mon = 5;
  else if(month == "JUL") temp_tm.tm_mon = 6;
  else if(month == "AUG") temp_tm.tm_mon = 7;
  else if(month == "SEP") temp_tm.tm_mon = 8;
  else if(month == "OCT") temp_tm.tm_mon = 9;
  else if(month == "NOV") temp_tm.tm_mon = 10;
  else if(month == "DEC") temp_tm.tm_mon = 11;

  // Process day:

  std::string day_str = t.substr(8, 2);
  
  if(day_str != "---")
    temp_tm.tm_mday = atoi(day_str.c_str());

  // Process seconds, minutes, hours, years:

  temp_tm.tm_hour = atoi(t.substr(11, 2).c_str());
  temp_tm.tm_min  = atoi(t.substr(14, 2).c_str());
  temp_tm.tm_sec  = atoi(t.substr(17, 2).c_str());
  temp_tm.tm_year = atoi(t.substr(20, 4).c_str()) - 1900;

  new_time = mktime(&temp_tm);

  return 0;
}

//=========================================================
//
//  indentString(...)
//
//=========================================================
inline std::string 
StringUtils::indentString(const std::string & src, size_t depth, bool skip_first_line)
{
  std::ostringstream s;
  std::ostringstream filler;
    
  if(depth)
    filler << std::setw(depth) << std::setfill(' ') << "";

  if(!skip_first_line)
    s << filler.str();

  for(size_t i = 0; i < src.length(); ++i)
  {
    s << src[i];

    if(src[i] == '\n')
      s << filler.str();
  }

  return s.str();
}

//=============================================================
//
//  lineWrapString(...)
//
//=============================================================
inline std::string 
StringUtils::lineWrapString(const std::string & src, size_t width)
{
  // Make sure width >= 1:
  if(width < 1)
    width = 1;

  // Set up the target stream and tokenize the input:

  std::ostringstream       tgt;
  std::vector<std::string> tokens(0);
  
  UTIL::StringUtils::splitOnWhitespace(src, tokens);
  
  // Loop across all tokens:
  
  size_t cur_line_width = 0;

  for(size_t i = 0; i < tokens.size(); ++i)
  {
    if(tokens[i].length() > width) // If the current token is itself greater than allowable width:
    {
      // Stick on a space or a newline:
      if(cur_line_width > 0)
      {
        if(cur_line_width < width) { tgt << " ";       cur_line_width++;   }
        else                       { tgt << std::endl; cur_line_width = 0; }
      }
      
      // Process it character by character:
      for(size_t j = 0; j < tokens[i].length(); ++j)
      {
        // Stick on a newline if necessary:
        if(cur_line_width == width) { tgt << std::endl; cur_line_width = 0; }
        
        // Append the j'th character of the token:
        tgt << tokens[i][j];
        cur_line_width++;
      }
    }
    else // If the current token is itself <= the allowable width
    {
      // If the appended current token exceeds allowable width:
      if(cur_line_width + 1 + tokens[i].length() > width)
      {
        tgt << std::endl;
        cur_line_width = 0;
      }
        
      // Append the token and increment the cur_line_width:
      tgt << (cur_line_width ? " " : "") << tokens[i];
      cur_line_width += (cur_line_width ? 1 : 0) + tokens[i].length();
    }
  }
  
  // Return string:
  return tgt.str();
}

//===========================================
//
//  stripLeadingExpZero
//
//===========================================
inline std::string 
StringUtils::stripLeadingExpZero(const std::string & s)
{
  std::string              pattern = "(.+e[\\+-])0(\\d\\d)";
  std::vector<std::string> results = UTIL::RegexUtils::matchPattern(pattern, s);

  return results.size() == 3 ? results[1] + results[2] : s;    
}


} // End namespace UTIL
} // End namespace SAGE

#endif
