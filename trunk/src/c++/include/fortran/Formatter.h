#ifndef __FORMATTER_H
#define __FORMATTER_H
//
//  Formatters -- Tells format of i/o
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   0.1 gcw Initial Implementation         Aug  1 96
//             1.0 gcw Revised for newest LSF (1.7?)  Nov 18 96
//             1.1 gcw Renaming and redoing certain   Oct 30 97
//                     bits
//
//  Copyright (c) 1997  R.C. Elston
//


#include <utility>
#include <stack>
#include <iostream>
#include "LSF/Attr.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

// A Formatter is a tool used by the tokenizer that knows the format of the
// i/o stream.  The tokenizer makes a request to the formatter, and the
// formatter responds with the next action to be performed.  This can
// be read/writing, moving elsewhere in the current line buffer, flushing
// the line (and reading the next) or it could be a message or output.
//
// The Formatter also tells the Tokenizer what sort of data is is parsing for,
// whether it's reached the end of its knowlege, or there is an error.
//
// Note, currently there are several fortran specific parts to the
// Tokenizer_Actions, which are the method of message passing.  This is
// due to Fortran being the only currently available Formatter and also
// that fortran has several strange requirements that are difficult to do
// in a generic manner.

// Tokenizer_Action - A common data structure encapsulating tokenizer
//     commands.

struct Tokenizer_Action
{
  enum Taction { move, read, flush, output, msg };
  enum d_type  { None,      Integer,   Float,     String, 
                 EndFormat, EndOfFile, Error };  
  enum flags   { clear, BZflag };

  Taction Action;        // Action to perform
  int     start, end;    // start and end position to read
  int     decimal;       // number of decimal places.
  d_type  data_type;     // Data type to read
  int     flag;          // Special control flags
  string  action_string; // Delimiters/error strings, etc.

  Tokenizer_Action() : flag(0) {}
};

// BaseFormatter - Doesn't do much. ;)  Is used as a base class for all the
//     different parsers we might want to develop.  Currently only fortran
//     Format statement style parsing is supported.

class BaseFormatter : public LSFBase
{
public:
  typedef Tokenizer_Action action;


  virtual void   reset ()           {                  }
  virtual action get_action()       { return action(); }
  virtual int    get_pos()    const { return -1;       }

protected:

  BaseFormatter() { }

};

// The stack_elt class is used internally by the format parser.  It is an
// element of the stack of recursions in the format statement.  The boolean
// operators are provided only for STL compatibility.

class stack_elt 
{
public:
  int start, rpt, end; // starting position, repetitions, end conditions
 
  stack_elt() {}
  stack_elt(int s, int r, int e) : start(s), rpt(r), end(e) {}
  bool operator==(const stack_elt &) const { return false; }
  bool operator!=(const stack_elt &) const { return false; }
  bool operator< (const stack_elt &) const { return false; }
  bool operator> (const stack_elt &) const { return false; }

};

// Currently the FortranFormatter is not strict.  It allows fortran statements
//     without parenthesis, and a variety of other common NoNos in strict
//     format statements (it allows whitespace, etc) We can change this, or
//     add a flag to handle it, if desired.
class FortranFormatter : public BaseFormatter
{
public:
 
  FortranFormatter()  { set_format(""); }
  ~FortranFormatter() { };
  
  void set_format (const string& s);

  virtual void   reset ( );
  virtual action get_action();
  virtual int    get_pos() const { return fmt_pos; }

protected:

#ifdef OLD_STL
  typedef stack< list<stack_elt> > Stack;
#else
  typedef stack< stack_elt, list<stack_elt> > Stack;
#endif  

  void ignore_ws(char * n = " \t") 
  {
    while ( fmt_pos < fmt_length && strchr(n, fmt_str[fmt_pos]) ) fmt_pos++;
  }

  void FORTRAN_ERROR   ( action& a, const char* err );
  bool next_format_elt ( action& a );
  bool next_bound_elt  ( action& a );
  bool get_value       ( action& a );
  bool get_tab         ( action& a );
  bool get_BZflag      ( action& a );
  bool get_end_fmt     ( action& a );
  bool get_quote       ( action& a );
  void get_fmt_digit   ();
  int  get_number      ();
  bool popstack        ();
  char fmt_chr         () const {return(fmt_pos<fmt_length)? fmt_str[fmt_pos]:'\0';}

  string       fmt_str;
  unsigned int fmt_pos;
  unsigned int fmt_length;
  int          state;
  Stack        st;
  int          flag;
  int          input_pos;
};

#endif
