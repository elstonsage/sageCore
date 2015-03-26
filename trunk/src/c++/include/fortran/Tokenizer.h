#ifndef __TOKENIZER_H
#define __TOKENIZER_H
//
//  Tokenizer -- Tokenizes Column delimited text based on a Formatters's
//      commands.  Can be easily upgraded to handle delimited text.
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   0.1 gcw Initial Implementation         Aug  1 96
//             1.0 gcw Revised for newest LSF (1.7?)  Nov 18 96
//             1.1 gcw Revised with intent to add     Oct 30 97
//                     the new Parser framework.
//
//  Copyright (c) 1997  R.C. Elston
//


#include <string>
#include <iostream>
#include "LSF/Attr.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "fortran/Formatter.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

// Tokenizer - takes an istream (or char* string name) and reads the stream
//     line by line.  It uses a Formatter to tell it the structure of the
//     file and how to parse it into tokens.  It can currently only handle
//     files where the data comes in a specific ordering, as the parser
//     does not get data passed to it.

class Tokenizer : public LSFBase
{
public:

  typedef Tokenizer_Action              action;
  typedef pair<action::d_type, AttrVal> token;

  // Note:  These are really incorrect.  Since Tokenizers now allow
  // both the Formatter and the istream to be replaced, we should allow that
  // here, and deal with correcting it as new things are entered.
  // This will simplify certain areas of the code, as we don't have to
  // create dummy Formatters and streams. - GCW 9711.13
  Tokenizer(BaseFormatter* _p, istream& _i);
  Tokenizer(BaseFormatter* _p, char* n    );

  ~Tokenizer();

// Formatter functions

  BaseFormatter* set_formatter(BaseFormatter* b)
  { if(b)                    { b->grab(); if(p) p->release(); p = b; }
    if(i && !i->fail() && p) clear();
    return p;
  }
  BaseFormatter* get_formatter() const { return p; }
  
// Istream functions

  istream& set_stream(istream& is) { i = &is; buffer = ""; return *i; }
  istream& get_stream() const      { return *i;          }

// Testing functions

  bool eol() const;    // Are we at end of current line?
  bool eof() const;    // Are we at the end of the file?

// String functions

  void  read_line();   // Read the next line into the buffer
  token get_token();   // Get the next token from the file

  string get_line() const { return buffer; }
  
  size_t line_count() const { return count; }

  size_t set_line_count(size_t s) { return count = s; }

protected:

  typedef Tokenizer_Action::d_type d_type;

  token read_token(action&, token&); 
  void  convert_string(string&, const string&, const d_type, const int );
  void  convert_float(action&, token&);

  istream*       i;
  BaseFormatter* p;
  string         buffer;
  bool           _eof;
  
  size_t         count;
};

// OutputTokenizer - takes an ostream and, upon being given tokens, puts
//     them into a buffer, which is written to the ostream.  Works pretty
//     much as the Tokenizer does.

class OutputTokenizer : public LSFBase
{
public:

  typedef Tokenizer_Action              action;
  typedef pair<action::d_type, AttrVal> token;

  OutputTokenizer(BaseFormatter* _p, ostream& _i);
  OutputTokenizer(BaseFormatter* _p, char* n    );

  ~OutputTokenizer();

  void  write_line();
  token put_token(const token&);

protected:

  void write_token(action&, const token&);
  void asterik(string& s)
  { for (string::iterator i = s.begin(); i != s.end(); i++) *i = '*'; }

  ostream*       o;
  BaseFormatter* p;
  string         buffer;
};

#endif
