#ifndef PARSE_OPS_H
#define PARSE_OPS_H
//
//  PARSE OPERATIONS 1.0 -- Low level stream parse operations
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation        Apr 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//



#include <string>
#include <iostream>
#include <utility>
#include <time.h>
#include <stdio.h>

using namespace std;

typedef pair<string, string> spair;

// Read white space and check next character is in delim (doesn't take char)
// Note: Does not count lines, so be careful using this
bool scan_delimit(FILE *i, const char *ws = " \t",
                  const char *delim = "\r\n");

int kill_ws(FILE *i, const char *ws = " \t", const char *eol = "\r\n");


string  getString(FILE *i, const char *delim = "",
                              const char *ws    = " \t\n\r" );

string getQString(FILE *i, const char *delim = "",
                              const char    *ws = " \t",
                              const char   *eol = "\n\r",
                              int *lines = NULL );

int        getInt(FILE *i, const char *delim = "",
                              const char *ws    = " \t\n\r" );

bool scan_delimit(istream &i, const char *ws = " \t",
                  const char *delim = "\r\n");

int kill_ws(istream &i, const char *ws = " \t", const char *eol = "\r\n");


string  getString(istream &i, const char *delim = "",
                              const char *ws    = " \t\n\r" );

string getQString(istream &i, const char *delim = "",
                              const char    *ws = " \t",
                              const char   *eol = "\n\r",
                              int *lines = NULL );

int        getInt(istream &i, const char *delim = "",
                              const char *ws    = " \t\n\r" );

void  printQuotedChars(ostream &o, const string &s);
int   printQuoted(ostream &o, string& s, int _pos = 0, int ind = 0,
                  int _start = 0);
                                
spair parseVar(const string& var);

int strcasecmp(const string& str1, const string &str2);

string toLower( const string &s );
string toUpper( const string &s );

string doub2str (double d,      int width = 0, int pres = -1, long flags = 0,
                 char fill = '\0');
string long2str (long d,        int width = 0, long flags = 0,
                 char fill = '\0');
string ptr2str  (const void* v, int width = 0, long flags = 0,
                 char fill = '\0');

string convert_time(time_t t, bool use_names = false);

double str2doub(const string &s);
long   str2long(const string &s);

float  normalize(float  val, float  cutoff = 1e-15);
double normalize(double val, double cutoff = 1e-15);


class AttrVal;

/// Strips whitespace from the given string.
/// \param s The string from which whitespace will be stripped
/// \param ws A string representing the set of characters to be consider 'whitespace'.
/// If unspecified, this function will strip the traditional notion of whitespace
/// (tabs, newlines, spaces).
/// \returns The stripped string
string strip_ws(const string& s, const char *ws = NULL);

AttrVal stringtol(const string&);
AttrVal stringtof(const string&);

// string_tokenizer class - provides const forward iteration through
// delimited strings.
//
// TODO: Add configurable quote characters, string assignment and 
//       other finishing touches

class string_tokenizer
{
public:
  class iterator;
  friend class iterator;

  typedef iterator const_iterator;
  string_tokenizer(const string& s = "", const string &d = ",", 
                                         const string &w =" \t\n\r")
    : st(s), delim(d), ws(w), skip_consecutive_delim(false),
      skip_leading_delim(false), skip_trailing_delim(false) { }

  ~string_tokenizer() { }
  
// Iterators  

  iterator begin() const { return iterator(*this, 0);           }
  iterator   end() const { return iterator(*this, st.size()+1); }

  const string &str()        const { return st;    }
  const string &whitespace() const { return ws;    }
  const string &delimiters() const { return delim; }

  void set_str(const string &s)        { st = s;    }
  void set_whitespace(const string &w) { ws = w;    }
  void set_delimiters(const string &d) { delim = d; }

  bool skip_consecutive_delimiters() const
  { 
    return skip_consecutive_delim;
  }
  bool skip_leading_delimiters() const
  { 
    return skip_leading_delim;
  }
  bool skip_trailing_delimiters() const
  { 
    return skip_trailing_delim;
  }
  void set_skip_consecutive_delimiters(bool skip = true) 
  { 
    skip_consecutive_delim = skip;
  }
  void set_skip_leading_delimiters(bool skip = true) 
  { 
    skip_leading_delim = skip;
  }
  void set_skip_trailing_delimiters(bool skip = true) 
  { 
    skip_trailing_delim = skip;
  }
  
protected:
  string st, delim, ws;
  bool skip_consecutive_delim;
  bool skip_leading_delim;
  bool skip_trailing_delim;

public:

  class iterator
  {
    //lint --e(1554,1555)
  public:
    iterator() : st(NULL), pos(0), first(false), last(false), last_delim('\0') { }
    iterator(const string_tokenizer&, size_t);
    iterator(const iterator& s)
       : st(s.st), pos(s.pos), first(s.first), last(s.last), 
         last_delim(s.last_delim), value(s.value)              { }
    
    ~iterator() { st = NULL; }
    
    iterator& operator= (const iterator& rhs)
    { 
      if(&rhs == this) return *this;
    
      st         = rhs.st; 
      value      = rhs.value; 
      first      = rhs.first;
      last       = rhs.last;
      pos        = rhs.pos;
      last_delim = rhs.last_delim;
      return *this; 
    }
    
    char last_delimiter() const { return last_delim; }

  // Dereference

    const string &operator*() const  { return value;  }
    const string *operator->() const { return &value; }

  // Comparisons

    bool operator==(const iterator& rhs) const
    { return st == rhs.st && pos == rhs.pos && first == rhs.first 
        && last == rhs.last; }

    bool operator!= (const iterator& rhs) const
    { return !(*this == rhs); }
    
  // iteration

    iterator& operator++();
    iterator  operator++(int) { iterator i(*this); ++(*this); return i; }

    iterator& operator--();
    iterator  operator--(int) { iterator i(*this); --(*this); return i; }
    
  protected:

    void find_end();
    void find_start(bool first_field = false);
    
    const string_tokenizer* st;
    size_t pos;
    bool first;
    bool last;
    char last_delim;
    string value;
  };
};

#endif
