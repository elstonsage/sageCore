#ifndef __DOTS_H
#define __DOTS_H

#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include <algorithm>
#include "LSF/parse_ops.h"

// Dot formatter base class

using namespace std;

/// The dot formatter is used for timing of long running processes.  In it,
/// you set a trigger count, which is the number of times it will be
/// triggered.  Several functions are provided in the virtual interface to
/// allow the formatter to produce messages or perform calculations as
/// needed.


class dot_formatter
{
public:

  inline dot_formatter();
  
  virtual inline ~dot_formatter() { }

  inline size_t set_trigger_count(size_t t);
  inline size_t get_trigger_count() const;
  
  // Returns true when done.  False otherwise.
  inline bool trigger();

  // Finish the current set of triggers
  inline bool finish();

  inline bool done() const;

  // If spammed, must reset everything.
  inline bool spammed();

protected:

  /// Do any prefix processing (on the first trigger) and display
  virtual bool prefix() = 0;

  /// Do any postfix processing (after last trigger) and display
  virtual bool postfix() = 0;

  /// Do any Resynchronization of display and current trigger values.
  /// Returns true if resync is successful, false if it fails.
  virtual bool resync() = 0;

  /// Inform the derived class that we're resetting to the beginning of the
  /// line.
  virtual bool reset() = 0;

  size_t my_trigger_count;
  size_t my_current_value;
};

// Text display Dot Formatter

class text_dot_formatter : public dot_formatter
{
public:

  typedef std::string string_type;

  inline text_dot_formatter(std::ostream&);

  // Set the number of characters output due to the trigger

  inline size_t set_output_character_count(size_t = 20);
  inline size_t get_output_character_count() const;

  inline char set_output_char(char = '.');
  inline char get_output_char() const;

  // Prefix and postfix setting

  // When the prefix string is modified, the current line is considered
  // Finished and a new line will be produced.  Prefix values may be
  // modified freely in between lines.
  inline const string_type& set_prefix(const string_type&);
  inline const string_type& get_prefix() const;

  inline size_t set_prefix_width(size_t t = 0, char = ' ');
  inline size_t get_prefix_width() const;
  inline char   get_prefix_char() const;

  inline const string_type& set_postfix(const string_type&);
  inline const string_type& get_postfix() const;
  
  // Set the width of the postfix and the fill character
  inline size_t set_postfix_width(size_t t = 0, char = ' ');
  inline size_t get_postfix_width() const;
  inline char   get_postfix_char() const;

protected:

  // Do any prefix processing and display
  virtual inline bool prefix();

  // Do any postfix processing and display
  virtual inline bool postfix();

  // Do any Resynchronization of display and current trigger values.
  virtual inline bool resync();

  virtual inline bool reset();

  std::ostream& my_stream;

  string_type my_prefix;
  string_type my_postfix;
  
  size_t my_prefix_width;
  size_t my_postfix_width;

  char my_prefix_char;
  char my_postfix_char;

  size_t my_output_count;
  
  char my_output_char;

  size_t my_current_width;
};

/// Displays a chart with time used and (approximate) times of completion.

class time_dot_formatter : public dot_formatter
{
public:

  typedef std::string string_type;

  inline time_dot_formatter(std::ostream&);

  inline size_t set_output_time_count(size_t = 6);
  inline size_t get_output_time_count() const;

  inline size_t set_max_time_interval(size_t = 3600);
  inline size_t get_max_time_interval() const;

  inline size_t set_min_time_interval(size_t = 5);
  inline size_t get_min_time_interval() const;

  inline size_t set_time_check_count(size_t = 100);
  inline size_t get_time_check_count() const;

  // Prefix and postfix setting

  // When the prefix string is modified, the current line is considered
  // Finished and a new line will be produced.  Prefix values may be
  // modified freely in between lines.
  inline const string_type& set_prefix(const string_type&);
  inline const string_type& get_prefix() const;

  inline size_t set_prefix_width(size_t t = 0, char = ' ');
  inline size_t get_prefix_width() const;
  inline char   get_prefix_char() const;

  inline const string_type& set_postfix(const string_type&);
  inline const string_type& get_postfix() const;
  
  // Set the width of the postfix and the fill character
  inline size_t set_postfix_width(size_t t = 0, char = ' ');
  inline size_t get_postfix_width() const;
  inline char   get_postfix_char() const;

protected:

  // Do any prefix processing and display
  virtual inline bool prefix();

  // Do any postfix processing and display
  virtual inline bool postfix();

  // Do any Resynchronization of display and current trigger values.
  virtual inline bool resync();

  virtual inline bool reset();

  inline void output_time(size_t elapsed, bool first);

  std::ostream& my_stream;

  size_t output_time_count;
  size_t max_time_interval;
  size_t min_time_interval;

  size_t time_check_count;

  string_type my_prefix;
  string_type my_postfix;
  
  size_t my_prefix_width;
  size_t my_postfix_width;

  char my_prefix_char;
  char my_postfix_char;

  time_t start_time;
  time_t current_time;
  
  size_t interval;
};

#include "util/dots.ipp"

#endif
