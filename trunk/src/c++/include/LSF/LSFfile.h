#ifndef __LSF_FILE_H
#define __LSF_FILE_H
//
//  LSF FILE OBJECTS 1.0 -- Low level logical structure building block objects
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation        Apr 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//


#include "LSF/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

class LSF_input : public LSFBase
{
public:
  typedef list<LSFBase *> path_t;
   
  LSF_input (istream &i,    ostream &o = cerr);
  LSF_input (const char *f, ostream &o = cerr);
  ~LSF_input();
  void input_to( LSFBase *o, bool typed = true );
  const path_t *path() { return &_path; }
  int lines() const    { return _lines; }
protected:
  typedef LSFDef::AList AList;
  void error(char *m=NULL);
  void comment();
  
  LSFBase *getLSF(bool typed = true);
  void  getPair( istream &i, const char *delim = "=,;#{}\n\r",
                             const char    *ws = " \t" );

  bool      free_file;
  istream  *lsf_file;
  ostream  *err;
  path_t   _path;
  AList    def;
  int      _lines;
};

// Load an LSF file and insert it into a new LSF base with name = desc
LSFBase *loadLSFfile(const string& file, const string& desc,
                     std::ostream& errors, bool typed = true);
                                     
class SimpleLSFDump : public LSFVisitor
{
public:
  SimpleLSFDump() {}
  void output(ostream &o, LSFBase *i);
  virtual void Process(LSFBase *c);

private:
  virtual void def(LSFBase *c);
  void output(LSFBase *c);
  void indent(int n);
  ostream *out;            
  int level;
};

class LSFDump : public LSFVisitor
{
public:
  LSFDump(char *name);
  LSFDump(ostream &o) { out = &o; free_stream = false; }
  virtual ~LSFDump() { if(free_stream && out) delete out; }
  
  void dump(LSFBase *);
  virtual void Process(LSFBase *c);

private:
  virtual void def(LSFBase *c);
  void output(LSFBase *c);
  void indent(int n);
  ostream *out;            
  int level;
  bool free_stream;
};

#endif
