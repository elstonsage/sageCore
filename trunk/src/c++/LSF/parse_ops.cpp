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
#include <algorithm>
#include <utility>
#include <iostream>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

#include <sstream>

#include "globals/config.h"
#include "LSF/Attr.h"
#include "LSF/parse_ops.h"

#ifndef _MAX
#  define _MAX(x, y) ( ( (x) < (y) )? (y) : (x) )
#endif

#define DEBUG(x) 

bool scan_delimit(FILE *i, const char *ws, const char *delim)
{
  if(!i || ferror(i) || feof(i)) return false;

  kill_ws(i, ws, "");

  if(ferror(i) || feof(i)) return false;

  unsigned int c = getc(i);  
  ungetc(c, i);
  return (strchr(delim,c) != NULL);
}

int kill_ws(FILE *i, const char *ws, const char *eol)
{
  if(!i || ferror(i) || feof(i)) return 0;

  int l = 0;
  unsigned int c;

//while(i.get(c))
  for(c=getc(i); !ferror(i) && !feof(i); c=getc(i))
  {
    if( strchr(ws, c) ) continue;
    else if( strchr(eol,c) ) { ++l; continue; }
    else break;
  }

  if(!feof(i))
    ungetc(c, i);

  return l;
}

// Get a string up to white space
string getString(FILE *i, const char *delim, const char *ws )
{
  const int bsize = 512;
  char buffer[512];
  int bcount = 0;
  string s;

  if(!i || feof(i) || ferror(i) )
    return s;

  char        c;
  const char* sel = NULL;

  for(c=getc(i); !ferror(i) && !feof(i); c=getc(i) )
  {      
    if( (sel=strchr(delim,  c)) != NULL ||
             strchr(ws,     c)  != NULL)
      break;

     if(bcount == bsize)
     {
       s.append(buffer, bcount);
       bcount = 0;
     }

     buffer[bcount++] = c;
  }

  if(bcount != 0)
    s.append(buffer, bcount);

  if(sel)
    ungetc(c,i);

  return s;
}

string getQString(FILE *i, const char *delim, const char    *ws,
                           const char   *eol, int *lines )
{
  char buffer[100];
  char *stop  = buffer;
  string s;

  if(!i || feof(i) || ferror(i) )
    return s;

  bool literal = false;
  bool quoted  = false;
  char        c;
  const char* sel = NULL;

  for(c=getc(i); !ferror(i) && !feof(i); c=getc(i) )
  {      
    sel = NULL;
    if( !quoted && !literal && ( (sel=strchr(delim,  c)) != NULL) )
      break;

    if(stop == buffer + 98)
    {
      *stop = '\0';
      s+=buffer;
      stop = buffer;
    }

    if(literal)                 
    { 
      literal = false;
      switch(c)
      {
        case 'n' : *stop++ = '\n'; break;
        case 't' : *stop++ = '\t'; break;
        default: if(!strchr(eol, c)) *stop++ = c; 
                 else if( lines ) (*lines)++;
                 break;
      }
      continue;
    }
    else if(c == '\\')            literal = true;
    else if(c == '\"' && quoted)  quoted  = false;
    else if(c == '\"')            quoted  = true;
    else if(strchr(eol, c))       break;
    else if(!quoted && strchr( ws, c )) /*skip white space*/ ;
    else *stop++ = c;
  }

  if(stop != buffer)
  {
    *stop = '\0';
    s+=buffer;
  }

  if(sel)
    ungetc(c,i);

  return s;
}

int getInt(FILE *i, const char *delim, const char *ws)
{
	return atoi( getString(i, delim, ws).c_str() );
}

bool scan_delimit(istream &i, const char *ws, const char *delim)
{
  if(!i || i.eof()) return false;

  kill_ws(i, ws, "");

  return (strchr(delim,i.peek()) != NULL);
}

int kill_ws(istream &i, const char *ws, const char *eol)
{
  if(!i || i.eof()) return 0;
  int l = 0;
  char c;

//while(i.get(c))
  for(c=i.get(); i; c=i.get())
  {
    if( strchr(ws, c) ) continue;
    else if( strchr(eol,c) ) { ++l; continue; }
    else break;
  }

  if(!i.eof())
    i.putback(c);

  return l;
}

// Get a string up to white space
string getString(istream &i, const char *delim, const char *ws )
{
  char buffer[100];
  char *stop  = buffer;
  string s = "";

  if(!i || i.eof())
    return s;

  char        c;
  const char* sel = NULL;
//while(i.get(c))
  for(c=i.get(); i; c=i.get() )
  {      
    if( (sel=strchr(delim,  c)) != NULL ||
             strchr(ws,     c)  != NULL)
      break;
     *stop++ = c;
     if(stop == buffer + 98)
     {
       *stop = '\0';
       s+=buffer;
       stop = buffer;
     }
  }

  if(stop != buffer)
  {
    *stop = '\0';
    s+=buffer;
  }

  if(sel)
    i.putback(c);

  return s;
}

string getQString(istream &i, const char *delim, const char    *ws,
                              const char   *eol, int *lines )
{
  char buffer[100];
  char *stop  = buffer;
  string s;

  if(!i || i.eof()) return s;

  bool literal = false;
  bool quoted  = false;
  char        c;
  const char* sel = NULL;

  for(c=i.get(); i; c=i.get())
  {      
    sel = NULL;
    if( !quoted && !literal && ( (sel=strchr(delim,  c)) != NULL) )
      break;

    if(stop == buffer + 98)
    {
      *stop = '\0';
      s+=buffer;
      stop = buffer;
    }

    if(literal)                 
    { 
      literal = false;
      switch(c)
      {
        case 'n' : *stop++ = '\n'; break;
        case 't' : *stop++ = '\t'; break;
        default: if(!strchr(eol, c)) *stop++ = c; 
                 else if( lines ) (*lines)++;
                 break;
      }
      continue;
    }
    else if(c == '\\')            literal = true;
    else if(c == '\"' && quoted)  quoted  = false;
    else if(c == '\"')            quoted  = true;
    else if(strchr(eol, c))       break;
    else if(!quoted && strchr( ws, c )) /*skip white space*/ ;
    else *stop++ = c;
  }

  if(stop != buffer)
  {
    *stop = '\0';
    s+=buffer;
  }

  if(sel)
    i.putback(c);
  return s;
}

int getInt(istream &i, const char *delim, const char *ws)
{
	return atoi( getString(i, delim, ws).c_str() );
}

void printQuotedChars(ostream &o, const string &s)
{
  for(unsigned int i=0; i < s.size(); ++i)
    switch( s[i] )
    {
      case '\r' : break;
      case '\n' : o << "\\n";  break;
      case '\t' : o << "\\t";  break;
      case '\\' : o << "\\\\"; break;
      case '\"' : o << "\\\""; break;
      default   : o <<  s[i];  break;
    }
}

int printQuoted(ostream &o, string& s, int pos, int mark, int m)
{
  // pos  - current position on the line.
  // mark - current tab amount
  // m    - position in the string to read from

  if(!o) return pos;

  int len = s.length();  // len - length of string
  int stop  = len;

  // While the whole string has not been written
  DEBUG(cout << "printQuoted(" << s << ")" << endl;)
  while (m < len)
  {
    if( m && pos != mark && len-m+pos > 75 && len-m > 1)
    {
      o << endl;
      for(pos=0; pos<mark; ++pos) o << ' ';
    }

    stop = len;

    // If the string is too long to fit on the current line
    if( pos + len + 2 - m > 76 )
    {
      stop = m+76-pos;
      if(stop > len)
        stop = len;
      for(int p = 5+m; p < len && pos + len - p + 2 < 77; ++p)
        if(strchr(" \t\n;.,*-+/|?!", s[p]))
          stop = p;
      if(stop <= m) stop = m+1;
    }

    o << '\"';                      // output \"
    printQuotedChars(o, s.substr(m, stop-m));
    o << '\"';
    pos += 2+stop-m;
    m = stop;

    if (m < len)
    {
      o << " \\";
      pos += 2;
    } 
  }
  return pos;
}

spair parseVar(const string& var)
{
  unsigned int s = var.length();

  unsigned int fs = 0;
  for (; (fs < s && !strchr(" \t[", var[fs])); fs++);
  if (fs == s) return spair();

  unsigned int bi = fs;
  for (; (bi < s && var[bi] != '[' && strchr(" \t", var[bi])); bi++);
  if (bi == s || var[bi] != '[') return spair();

  for (bi++; (bi < s && strchr(" \t", var[bi])); bi++);
  if (bi == s) return spair();

  unsigned int ei = bi + 1;
  for (; (ei < s && !strchr(" \t]", var[ei])); ei++);
  if (ei == s) return spair();

  unsigned int rb = ei;
  for (; (rb < s && var[rb] != ']' && strchr(" \t", var[rb])); rb++);
  if (rb != s-1 || var[rb] != ']') return spair();

  return spair(var.substr(0,fs), var.substr(bi, ei-bi) );
}

int strcasecmp (const string& str1, const string &str2)
{
  return strcasecmp( str1.c_str(), str2.c_str() );
}

string toLower( const string &s )
{
  string ls = s;
  for(unsigned int i=0; i<s.length(); i++) ls[i] = tolower( s[i] );
  return ls;
}

string toUpper( const string &s )
{
  string us = s;
  for(unsigned int i=0; i<s.length(); i++) us[i] = toupper( s[i] );
  return us;
}

string doub2str (double d, int width, int pres, long flags, char fch)
{
  if(width < 0) width = 0;

  bool b = (pres < 0);                  // pres < 0 indicates variable pres.

  if(pres < 0 && width) pres = width;

  while(true)
  {
    ostringstream o;

    if(flags)     o << setiosflags( (ios_base::fmtflags)flags );
    if(width)     o << setw(width);
    if(fch)       o << setfill(fch);
    if(pres >= 0) o << setprecision( (unsigned) pres);

    o << d;
    int size = o.str().size();

    if(b && width &&  size > width)
    {
      if(pres) pres -= (size - width); // pres must exist and get smaller
      else
	  {
		  return o.str();
    
	  }
    }
	else 
    {
	  return o.str();
	}
    
	if(pres < 0) pres = 0;  // next time abort
  }
}

string long2str (long l, int width, long flags,
                 char fch)
{

  ostringstream o;

  if(flags) o << setiosflags( (ios_base::fmtflags)flags );
  if(width) o << setw(width);
  if(fch)   o << setfill(fch);

  o << l;

  return o.str();
}

string ptr2str (const void* v, int width, long flags, char fch)
{
  ostringstream o;

  if(flags) o << setiosflags( (ios_base::fmtflags)flags);
  if(width) o << setw(width);
  if(fch)   o << setfill(fch);

  o << v;

  return o.str();
}

string convert_time(time_t t, bool use_name)
{
  long seconds = t % 60;
  long minutes = t / 60 % 60;
  long hours   = t / 3600 % 24;
  long days    = t / 3600 / 24;

  string my_time;

  int first = 0;

  if(days)
  {
    my_time += long2str(days);

    if(use_name) 
      my_time += 'd';
    else
      my_time += ':';
    first = 1;
  }

  if(hours || first)
  {
    my_time += long2str(hours, 2 * first, 0, '0');

    if(use_name) 
      my_time += 'h';
    else
      my_time += ':';

    if(!first) first = 2;
  }

  if(first != 1 || !use_name)
  {
    my_time += long2str(minutes, 2 * (bool) first, 0, '0');

    if(use_name) 
      my_time += 'm';
    else
      my_time += ':';

    if(!first) first = 3;
  }

  if(first == 0 || first == 3 || !use_name)
  {
    my_time += long2str(seconds, 2 * (bool) first, 0, '0');

    if(use_name) 
      my_time += 's';
  }

  return my_time;
}

string strip_ws(const string& s, const char *ws)
{
  size_t b, e;

  if(s.length() == 0) return string();

  if(ws)
  {
    b = s.find_first_not_of(ws);
    e = s.find_last_not_of(ws);
  }
  else
  {
    for(b=0; b < s.length() && isspace(s[b]); ++b);
    for(e=s.length()-1; isspace(s[e]); --e);
  }

  if(e < b || b == (size_t)-1) return string();

  return s.substr(b,e-b+1);
}

AttrVal stringtol(const string& val)
{
  char *end_ptr = NULL;
  char *c = new char[val.length()+1];
  strcpy( c, val.c_str() );
  c[val.length()] = '\0';
 
  AttrVal v = (int) strtol(c, &end_ptr, 0);
  if(end_ptr && end_ptr != &c[val.length()])   // Invalid Input
  {
    delete[] c;
    return AttrVal();
  }
  delete[] c;
  return v;
}

AttrVal stringtof(const string& val)
{
  char *end_ptr = NULL;
  const char *begin = val.c_str();
 
  double v = strtod(begin, &end_ptr);
  if(end_ptr == begin)   // Invalid Input -- cannot create anything
    return AttrVal();

  return AttrVal(v);
}

double str2doub(const string& val)
{
  char *end_ptr = NULL;
  const char *begin = val.c_str();
 
  double v = strtod(begin, &end_ptr);
  if(end_ptr == begin)   // Invalid Input -- cannot create anything
    return numeric_limits<double>::quiet_NaN();
  return v;
}

long str2long(const string& val)
{
  char *end_ptr = NULL;
  const char *begin = val.c_str();
 
  long v = strtol(begin, &end_ptr, 10);
  if(end_ptr == begin)   // Invalid Input -- cannot create anything
    return 0;
  return v;
}

float normalize(float f, float cutoff)
{
  if(fabs(f) < cutoff) return 0.;

  return f;
}

double normalize(double d, double cutoff)
{
  if(fabs(d) < cutoff) return 0.;

  return d;
}

// String Tokenizer stuff

string_tokenizer::iterator::iterator(const string_tokenizer& s, size_t i)
{
  st = &s;
  pos = i;

  value.resize(0);
  first = last = false;
  last_delim = '\0';

  find_start( true );
  find_end();
}

string_tokenizer::iterator& string_tokenizer::iterator::operator++()
{
  value.resize(0);
  last_delim = '\0';

  find_start();
  find_end();
  return *this;
}

void string_tokenizer::iterator::find_start(bool first_field)
{
  const string &s=st->st;

  // Catch out of bounds state
  if(pos >= s.size())
    return;

  // Construct shorthand notation for the following logic
  bool l = st->skip_leading_delimiters();
  bool t = st->skip_trailing_delimiters();
  bool c = st->skip_consecutive_delimiters();

  // See state table below for justification
  // Intuitively, no delimiter eliding is required so we return.
  if( !l && !t && !c )
    return;

  // Scan string for last consecutive delimiter
  size_t i = pos;
  for(; i < s.size() && strchr(st->delim.c_str(), s[i]) != NULL; ++i);

  // If no characters after the current are delimiters, then return
  if(i==pos) return;

  bool f = first_field;
  bool e = (i == s.size());

  // The next section of code implements the rules for when delimiters are
  // elided based on a set of state variables.  The logic looks fairly
  // obscure, but is based on the following state table and the previously
  // defined 5 variables (f,e,l,t,c).  Some cases are non-trivial to derive
  // since they rely on related behaviors of the 'first' and 'last'
  // variables.
  //
  //        \  L| 0 0 0 0 1 1 1 1
  //         \ T| 0 0 1 1 1 1 0 0
  //        FE\C| 0 1 1 0 0 1 1 0
  //        ----|-----------------    Coverage of R:
  //        0 0 | R     R R     R       1) !T !C !F
  //        0 1 | R             R       2)  T !C !F !E
  //        1 1 | R R R R     R R       3) !L  F
  //        1 0 | R R R R               4)  F  E !L !T
  //        
  //        R     = return without action
  //        blank = skip consecutive delimiters         

  if( (!t && !c && !f) || (t && !c && !f && !e) || (!l && f) || (f && e && !l && !t))
    return;

  // Update the current position if in bounds
  if(i <= s.size())
    pos=i;

  // Invoke the trailing field flag if we do not elide trailing fields and
  // a delimiter character is at the end of the string.
  if( !t && e )
  {
    last = true;
    pos=i-1;
  }
}

void string_tokenizer::iterator::find_end() 
{
  const string &s = st->st;

  last_delim = '\0';

  if(pos >= s.size())
  {
    pos = s.size();
    first = last = true;
    return;
  }

  bool literal    = false;
  bool quoted     = false;
  bool last_quote = false;
  char c;

  const int bsize = 512;
  char buffer[512];
  int bcount = 0;
  size_t size = s.size();  

  while(pos < size)
  {
    c=s[pos++];

    if(    !quoted 
        && !(c == '"' && last_quote) 
        && !literal 
        && strchr(st->delim.c_str(),  c) != NULL )
    {
      last_delim = c;

      if(last)
        last_delim = '\0';

      // Invoke last field rule (last character is a delimiter)
      if(!last && pos == s.size())
      {
        --pos;
        last = true;
      }

      break;
    }

    if(bcount == bsize)
    {
      value.append(buffer, bcount);
      bcount = 0;
    }

    if(literal)                 
    { 
      literal = false;
      switch(c)
      {
        case 'n' : buffer[bcount++] = '\n'; break;
        case 't' : buffer[bcount++] = '\t'; break;
        default  : buffer[bcount++] = c;    break;
      }
      continue;
    }
    else if(c == '\\')
    {
      last_quote       = false;
      literal          = true;
    }
    else if(c == '"')
    {
      if(last_quote)
      { 
        last_quote = false; 
        quoted     = !quoted; 
        buffer[bcount++] = c; 
      }
      else if(quoted)
      {
        quoted     = false;
        last_quote = true;
      }
      else
      {
        // last_quote = true; // the first quote is never an internal quote
        quoted     = true;
      }
    }
    else 
    {
      last_quote       = false;
      buffer[bcount++] = c;
    }
  }

  if(literal)
    last_delim = '\\';

#if REPORT_ERROR
  if(quoted)
  // report error in file format
#endif

  if(bcount != 0)
    value.append(buffer, bcount);

  if(st->ws.size())
    value = strip_ws(value, st->ws.c_str());
  else
    value = strip_ws(value);
}

#undef DEBUG
