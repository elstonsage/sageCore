#ifndef __ERROR_HANDLER_H
#define __ERROR_HANDLER_H
//
//  LSF Error Handler 0.1 -- Basic Error handling until we can replace it
//		with somehting better.
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu
//
//  History:   0.01 gcw  Initial implementation               Aug 27 1996
//       
//  Copyright (c) 1996  R.C. Elston
//


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "LSF/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

namespace SAGE
{

enum error_t { Warning, Information, Error };

class ErrorHandler;

extern ErrorHandler* GlobalError;

// The error handler maintains a vector of ostreams and directs error messages
//   to the appropriate ostreams.  Since ostreams can share the same value,
//   it is possible to direct an error to more than one ostream. Value 1
//   defaults to standard out, value 2 is standard error.  Commonly, the
//   'summary file' will be value 3, and the output file is value 4. 
//   Summary and Output correspond to fort.21 and fort.22 respectively.

struct ErrorStruct
{
public:
	ErrorStruct() : o(NULL) {}
	ErrorStruct(ostream* _o, int _v) : v(_v)
	{ if (_o && _o->good()) o = _o; else o = NULL; }

	ostream* o;
	int      v;

	bool operator<(const ErrorStruct &) const { return false; }
    bool operator==(const ErrorStruct &) const { return false; }
    bool operator!=(const ErrorStruct &) const { return true; }
};


class ErrorHandler : public LSFBase
{
public:

  ErrorHandler(const string& s) : pname(s) 
  { insert(&cerr, 2); insert(&cout, 1); }

  ~ErrorHandler() {}
  
  void insert(ostream* o, int i) { estreams.push_back(ErrorStruct(o, i)); }

#if 0
  void erase(int i) 
  {
    for (size_t j = 0; j < estreams.size(); j++)
      if (estreams[j].v == i) { estreams.erase(&estreams[j]); j--; }
  }
#endif
  
  void error(int i, error_t e, const string& s)
  {
    for(size_t j = 0; j < estreams.size(); j++)
      if(estreams[j].v == i) 
        out_error(estreams[j].o, e, s);
  }

  ostream* out(int i, const string& s = "")
  {
    ostream* o = NULL;
    for(size_t j = 0; j < estreams.size(); j++)
      if(estreams[j].v == i)
      { 
        if(!o)       o = estreams[j].o;
        if(s.size()) (*estreams[j].o) << s;
      }
    return o;
  }

protected:
  typedef vector<ErrorStruct>        err_vector;

  void out_error(ostream* o, error_t e, const string& s);

  err_vector estreams;
  string     pname;
};

inline void ErrorInit(const string& s) { GlobalError = new ErrorHandler(s); }

}

#endif
