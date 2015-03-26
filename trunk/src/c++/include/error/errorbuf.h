#ifndef ERROR_BUF_H
#define ERROR_BUF_H

//==================================================================
//  File:       errorbuf.h
//
//  Purpose:    This header file defines stream buffer types for
//              error streams derived from basic_errorstream.
//
//  Author:     Kevin Jacobs
//
//  History:    Version 0.10
//
//  Copyright (c) 1998 R.C. Elston.  All Rights Reserved
//==================================================================

#include <iostream>
#include <sstream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <string>
#include <time.h>
#include "globals/config.h"
#include "boost/smart_ptr.hpp"

#define DEBUG(x)

namespace SAGE {

//==================================================================
//=                      Forward Declarations                      =
//==================================================================

template<class charT, class traits>
class basic_errorstream;


//==================================================================
//=            Definition of error priority levels                 =
//==================================================================

//lint ++flb We put this into library mode for the definition of the enummerations
//           This is to remove messages about unused enumerated values.
//           We may want to eventually make changes to our lint configuration to
//           understand our modules.

enum error_priority
{
  no_priority=-1,
  debug=1,
  information=50,
  notice=60,
  warning=70,
  error=80,
  critical=90, 
  fatal=100
};

enum error_restriction
{
  r_none,
  r_eq,
  r_ne,
  r_lt,
  r_gt,
  r_le,
  r_ge,
  r_range
};

//lint --flb Turn off library mode.

//==================================================================
//             Definition of basic error buffer
//------------------------------------------------------------------
//  Defines a base class for error string buffers.
//==================================================================

template<class charT,class traits=std::char_traits<charT> >
class basic_errorbuf : public std::basic_stringbuf<charT,traits>

{
  public:
    //lint -e{1516} Must redefine this, since it's private in the orginal
    typedef std::basic_string<charT,traits> string_type;

    typedef error_priority error_priority_type;
    //typedef basic_errorbuf<charT,traits> buf_type;

    struct error_info;

    explicit basic_errorbuf();
    basic_errorbuf(const basic_errorbuf &sb);
    virtual ~basic_errorbuf();

//    void str(const string_type& str);
//    string_type str() const;

    error_info &get_info();

    void prefix(const string_type &pre);
    string_type prefix() const;
    void suffix(const string_type &pre);
    string_type suffix() const;

    void location(const string_type &file, size_t line = ((size_t)-1));
    void filename(const string_type &file);
    string_type filename() const;
    void linenumber(size_t l);
    size_t linenumber() const;

    void priority(error_priority_type p);
    error_priority_type priority() const;

    void offset(size_t p);
    size_t offset() const;

    void indent(size_t p);
    size_t indent() const;

    void line_width(size_t p);
    size_t line_width() const;

    void set_raw_mode();
    bool raw_mode() const;

    void set_cooked_mode();
    bool cooked_mode() const;

    struct error_info
    {
        friend class basic_errorbuf<charT,traits>;

        size_t _offset;
        size_t _indent;
        size_t _width;
        size_t _line;
        bool   _raw;

        string_type file;
        error_priority_type prior;
        string_type pre, suf;

        void assign(const error_info &e)
	{
          _offset = e._offset;
          _indent = e._indent;
          _width  = e._width;
          _line   = e._line;
          _raw    = e._raw;
          prior   = e.prior;
          file    = e.file;
          pre     = e.pre;
          suf     = e.suf;
	};

      private:
        //lint -e{1704} Private access desiered
        error_info()
        {
          prior   = no_priority;
          _width  = 79;
          _offset = 0;
          _indent =((size_t)-1);
          _line   =((size_t)-1);
          _raw    = false;
        }

        //lint -e{1704} Private access desiered
        error_info(const error_info &e)
        {
          assign(e);
        }
        error_info &operator=(const error_info &rhs);

    };

  protected:
    virtual int sync();
    error_info info;
};

//=================================================================
//           Implementation of basic error buffer class
//=================================================================

/******************************************************************
 *       Implementation of basic error buffer (non-INLINE)        *
 ******************************************************************/

template<class charT,class traits>
basic_errorbuf<charT,traits>::
  basic_errorbuf() :
    std::basic_stringbuf<charT,traits>(std::ios_base::out|std::ios_base::app)
{
  DEBUG(std::cout << "new basic_errorbuf() at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errorbuf<charT,traits>::basic_errorbuf
(const basic_errorbuf<charT,traits> &sb) :
   std::basic_stringbuf<charT,traits>(sb)
{
  str(sb.str());

  info.prior   = sb.info.prior;
  info.suf     = sb.info.suf;
  info.pre     = sb.info.pre;
  info._offset = sb.info._offset;
  info._indent = sb.info._indent;
  info._width  = sb.info._width;
  info._raw    = sb.info._raw;
}

template<class charT,class traits>
basic_errorbuf<charT,traits>::
  ~basic_errorbuf()
{
  DEBUG(std::cout << "basic_errorbuf being synced before delete at " << (void*) this << std::endl;)
  //lint -e{534}  We also ignore it's return, because we don't care.
  basic_errorbuf<charT,traits>::sync();
  DEBUG(std::cout << "basic_errorbuf being deleted at " << (void*) this << std::endl;)
}

template<class charT,class traits>
int basic_errorbuf<charT,traits>::sync()
{
  return -1;
}

     /*-----------------------------------------------------*
      *    Implementation of basic error buffer (inline)    *
      *-----------------------------------------------------*/
/*
template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::str(const string_type &s)
{
  std::basic_stringbuf<charT,traits>::str(s);
}


template<class charT,class traits>
inline basic_errorbuf<charT,traits>::string_type
basic_errorbuf<charT,traits>::str() const
{ return std::basic_stringbuf<charT,traits>::str(); }
*/
template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::prefix(const string_type &p)
{ info.pre = p; }

template<class charT,class traits>
inline typename basic_errorbuf<charT,traits>::string_type
basic_errorbuf<charT,traits>::prefix() const
{ return info.pre; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::suffix(const string_type &s)
{ info.suf = s; }

template<class charT,class traits>
inline typename basic_errorbuf<charT,traits>::string_type
basic_errorbuf<charT,traits>::suffix() const
{ return info.suf; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::location(const string_type &file,
                                                   size_t line)
{
  info.file  = file;
  info._line = line;
}

template<class charT,class traits>
inline typename basic_errorbuf<charT,traits>::string_type
basic_errorbuf<charT,traits>::filename() const
{ return info.file; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::filename(const string_type &file)
{ info.file = file; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::linenumber(size_t l)
{ info._line = l; }

template<class charT,class traits>
inline size_t basic_errorbuf<charT,traits>::linenumber() const
{ return info._line; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::priority(error_priority_type p)
{ info.prior = p; }

template<class charT,class traits>
inline typename basic_errorbuf<charT,traits>::error_priority_type
basic_errorbuf<charT,traits>::priority() const
{ return info.prior; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::indent(size_t i)
{ info._indent = i; }

template<class charT,class traits>
inline size_t basic_errorbuf<charT,traits>::indent() const
{ return info._indent; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::offset(size_t i)
{ info._offset = i; }

template<class charT,class traits>
inline size_t basic_errorbuf<charT,traits>::offset() const
{ return info._offset; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::line_width(size_t i)
{ info._width = i; }

template<class charT,class traits>
inline size_t basic_errorbuf<charT,traits>::line_width() const
{ return info._width; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::set_raw_mode()
{ info._raw = true; }

template<class charT,class traits>
inline bool basic_errorbuf<charT,traits>::raw_mode() const
{ return info._raw; }

template<class charT,class traits>
inline void basic_errorbuf<charT,traits>::set_cooked_mode()
{ info._raw = false; }

template<class charT,class traits>
inline bool basic_errorbuf<charT,traits>::cooked_mode() const
{ return !info._raw; }

template<class charT,class traits>
inline typename basic_errorbuf<charT,traits>::error_info &
basic_errorbuf<charT,traits>::get_info()
{
  //lint -e{1536} Exposing info this way is not good, but leave for now.
  return info;
}


//==================================================================
//              Definition of basic error stream buffer
//-----------------------------------------------------------------
//  Defines an error buffer class for pushing errors to an ostream
//==================================================================

template<class charT,class traits=std::char_traits<char> >
class basic_errorstreambuf : public basic_errorbuf<charT,traits>
{
  //lint --e{1712} No default constructor ok
  
  public:
    typedef basic_errorbuf<charT,traits> base_type;

    typedef typename base_type::string_type         string_type;
    typedef typename base_type::error_priority_type error_priority_type;

    //typedef basic_errorstreambuf<charT,traits> buf_type;

    explicit basic_errorstreambuf(std::ostream &o);
    explicit basic_errorstreambuf(std::ostream &o, const string_type& s);
    basic_errorstreambuf(const basic_errorstreambuf &sb);

    virtual ~basic_errorstreambuf();

    std::ostream &out_stream() const ;

  protected:
    virtual int sync();
    virtual void push_error();
    virtual size_t print_line_wrap(const string_type &line,
                                   size_t indent, size_t width);
    virtual string_type expand_string(const string_type &s) const;
    virtual string_type expand_option(char o) const;
    std::ostream &out;
};

//=================================================================
//     Implementation of basic error stream buffer (INLINE)
//=================================================================

template<class charT,class traits>
inline std::ostream &
basic_errorstreambuf<charT,traits>::out_stream() const
{
  //lint -e{1536} Exposing out is sometimes needed.
  return out;
}

/******************************************************************
 *   Implementation of basic error stream buffer (non-INLINE)     *
 ******************************************************************/

template<class charT,class traits>
basic_errorstreambuf<charT,traits>::basic_errorstreambuf(std::ostream &o) :
        basic_errorbuf<charT,traits>(), out(o)
{
    DEBUG(std::cout << "new basic_errorstreambuf() at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errorstreambuf<charT,traits>::
basic_errorstreambuf(std::ostream &o, const string_type& s) :
        basic_errorbuf<charT,traits>(), out(o)
{
  DEBUG(std::cout << "new basic_errorstreambuf() at " << (void*) this << std::endl;)
  str(s);
}

template<class charT,class traits>
basic_errorstreambuf<charT,traits>::
basic_errorstreambuf(const basic_errorstreambuf<charT,traits> &sb) :
       basic_errorbuf<charT,traits>(sb), out(sb.out_stream())
{
  str(sb.str());
}

template<class charT,class traits>
basic_errorstreambuf<charT,traits>::~basic_errorstreambuf()
{
  DEBUG(std::cout << "basic_errorstreambuf being synced before delete at " << (void*) this << std::endl;)
  //lint -e{534}  We also ignore it's return, because we don't care.
  basic_errorstreambuf<charT,traits>::sync();
  DEBUG(std::cout << "basic_errorstreambuf being deleted at " << (void*) this << std::endl;)
}

template<class charT,class traits>
typename basic_errorstreambuf<charT,traits>::string_type
basic_errorstreambuf<charT,traits>::expand_option(char o) const
{
  string_type opt;

  switch( o )
  {
    case '%':
      opt += '%';
      break;
    case 'D':
    case 'd':
      {
        time_t now = time(NULL);
        char *now_str = asctime( localtime( &now ) );
        size_t j = strlen(now_str) - 1;
        for(; (int)j >= 0 && isspace(now_str[j]) ; ++j)
          now_str[j] = '\0';
        opt += now_str;
      }
    break;
    case 'T':
    case 't':
      {
        time_t now = time(NULL);
        struct tm *now_tm = localtime( &now );
        char buf[10];
        sprintf(buf, "%02d:%02d:%02d", now_tm->tm_hour, now_tm->tm_min,
                                       now_tm->tm_sec);
        opt += buf;
      }
      break;
    case 'p':
      {
        int p = basic_errorbuf<charT,traits>::priority();
        if     (p < debug       ) break;
        else if(p < information ) opt += "D";
        else if(p < notice      ) opt += "I";
        else if(p < warning     ) opt += "N";
        else if(p < error       ) opt += "W";
        else if(p < critical    ) opt += "E";
        else if(p < fatal       ) opt += "C";
        else                      opt += "F";
      }
      break;
    case 'P':
      {
        int p = basic_errorbuf<charT,traits>::priority();
        if     (p < debug       ) opt += "None";
        else if(p < information ) opt += "Debug";
        else if(p < notice      ) opt += "Info";
        else if(p < warning     ) opt += "Notice";
        else if(p < error       ) opt += "Warning";
        else if(p < critical    ) opt += "Error";
        else if(p < fatal       ) opt += "Critical";
        else                      opt += "Fatal";
      }
      break;
    case 'n':
      {
        std::ostringstream oss;
        oss << (int)basic_errorbuf<charT,traits>::priority();
        opt += oss.str();
      }
      break;
    case 'w':
      {
        std::ostringstream oss;
        oss << basic_errorbuf<charT,traits>::line_width();
        opt += oss.str();
      }
      break;
    case 'F':
      if(!basic_errorbuf<charT,traits>::filename().size())
        opt += "none";
      break;
    case 'f':
      opt += basic_errorbuf<charT,traits>::filename();
      break;
    case 'L':
      if(basic_errorbuf<charT,traits>::linenumber() == (size_t)-1)
        opt += "none";
      break;
    case 'l':
      if(basic_errorbuf<charT,traits>::linenumber() != (size_t)-1)
      {
        std::ostringstream oss;
        oss << basic_errorbuf<charT,traits>::linenumber();
        opt += oss.str();
      }
      break;
    default:
      opt += '%';
      opt += o;
      break;
  }
  return opt;
}

template<class charT,class traits>
typename basic_errorstreambuf<charT,traits>::string_type
basic_errorstreambuf<charT,traits>::expand_string(const string_type &s) const
{
  if(!s.size()) return string_type();

  // Search strings for % characters to expand
  string_type ren;
  for(size_t i=0; i < s.size(); ++i)
  {
    switch( s[i] )
    {
      case '%':
      {
        if( ++i >= s.size() )
        {
          ren += '%';
          break;
        }
        ren += expand_option( s[i] );
        break;
      }

      default:
        ren += s[i];
        break;
    }
  }
  return ren;
}

template<class charT,class traits>
void basic_errorstreambuf<charT,traits>::push_error()
{
  DEBUG(std::cout << "basic_errorstreambuf::push_error() at " << (void*) this << std::endl;)
  if(!basic_errorbuf<charT,traits>::str().length()) return;
  if (!out || !out.good()) return;

  if( basic_errorbuf<charT,traits>::raw_mode() )
  {
    out << basic_errorbuf<charT,traits>::str() << std::flush;
    return;
  }

  size_t pos = 0;
  size_t ind = basic_errorbuf<charT,traits>::indent();
  string_type s = expand_string(basic_errorbuf<charT,traits>::str());

  // Expand and print prefix if it exists
  string_type fix = expand_string(basic_errorbuf<charT,traits>::prefix());
  if(fix.size())
  {
    if((int)ind<0) ind = fix.size();
    else ind += fix.size();

    // Print initial offset, prefix and update the line position
    for (size_t j = 0; j < basic_errorbuf<charT,traits>::offset(); ++j) out << ' ';
    
    out << fix;
    pos = basic_errorbuf<charT,traits>::offset() + ind;
  }

  // Compute the working indentation (based on prefix size and offset)
  if((int)ind<0)
    ind = 0;

  ind += basic_errorbuf<charT,traits>::offset();

  size_t l_stop;
  size_t l_start;
  size_t l_width = basic_errorbuf<charT,traits>::line_width();

  // Handle pathalogical case where line width-indent <= 0
  if(l_width > 0)
  {
    if(l_width <= ind)
      l_width = basic_errorbuf<charT,traits>::line_width() / 2 + 4;
    else
      l_width -= ind;
  }

  // Print each chunk of the string delimited by line break characters
  for (l_start = 0; l_start < s.length(); l_start = l_stop + 1)
  {
    // Find the beginning of the next chunk skipping whitespace
    for( ; l_start < s.length() && isspace(s[l_start]); ++l_start)
      //lint -e{722} Do nothing
      ;

    // Find line break characters to end the chunk
    l_stop  = s.find_first_of("\r\n",l_start);
    if((int)l_stop < 0)
      l_stop = s.length();

    size_t l = l_stop-l_start;
    if( l > 0 )
    {
      // Update line position and indent
      if( pos != 0 && pos != ind )
        out << '\n';
      if( pos != ind )
      {
        for (size_t j = 0; j < ind; ++j) out << ' ';
        pos = ind;
      }

      // Print line with special routine if it requires word wrap
      if(l_width > 0 && l >= l_width)
        pos = print_line_wrap( s.substr(l_start, l), ind, l_width );
      else
      {
        out << s.substr(l_start, l);
        pos = ind + l;
      }
    }
  }

  // Expand and print suffix
  fix = expand_string(basic_errorbuf<charT,traits>::suffix());
  if(fix.size())
  {
    // If suffix is too long then wrap it to the next line
    if( fix.size() + pos >= basic_errorbuf<charT,traits>::line_width() )
    {
      out << '\n';
      for (size_t j = 0; j < ind; ++j) out << ' ';
      pos = ind;
    }
    out << fix;
    pos += fix.size();
  }

  // Terminate with a new line if necessary and flush buffer
  if( pos != 0 )
    out << '\n';
  out << std::flush;
}

template<class charT,class traits>
size_t basic_errorstreambuf<charT,traits>::
print_line_wrap(const string_type &line, size_t ind, size_t width)
{
  if(!line.length()) return ind;

  size_t pos = 0;
  size_t l_break;
  size_t l_pos = (size_t)-1;

  // This should never happen
  if(ind >= width)
  {
    out << line;
    l_pos = ind + line.length();
  }

  while (pos < line.length())
  {
    // Skip leading white space
    if( isspace(line[pos]) )
    {
      ++pos;
      continue;
    }

    // Print indent for every line after the first
    if(l_pos != (size_t)-1)
    {
      out << out.widen('\n');
      for (size_t j = 0; j < ind; ++j) out << ' ';
      l_pos = ind;
    }

    // Compute upper bound of the current position that can be printed
    l_break = pos + width;

    // Check to see if the string is really that long
    if (l_break < line.length())
    {
      // If so then find a non-alphanumeric to break the string after
      size_t l_min = pos + std::min( (unsigned long) 2, (unsigned long)width/2 );

      size_t space;
      size_t punct;
      for( space = l_break; space > l_min && !isspace(line[space]); --space)
        //lint -e{722} Do nothing
        ;
      for( punct = l_break; punct > space &&  isalnum(line[punct]); --punct)
        //lint -e{722} Do nothing
        ;

      if( space > l_min && (size_t)abs( (int)space - (int)punct) < (width/10)+1 )
        l_break = space;
      else if( punct > l_min )
      {
        if( punct < l_break && strchr(".,!?)}]", line[punct] ) )
          l_break = punct + 1;
        else
          l_break = punct;
      }
    }

    if( l_break - pos > 0 )
    {
      // Output the current chunk
      out << line.substr(pos, l_break-pos);

      // Record the current line position
      l_pos = ind + l_break-pos;
    }

    // Find the beginning of the next chunk
    for( pos = l_break; pos < line.length() && isspace(line[pos]); ++pos)
      //lint -e{722} Do nothing
      ;
  }
  return l_pos;
}

template<class charT,class traits>
int basic_errorstreambuf<charT,traits>::sync()
{
  DEBUG(std::cout << "basic_errorstreambuf::sync() at " << (void*) this << std::endl;)
  if(!out)
    return -1;

  push_error();

  // Do whatever other stringbufs do to sync
  if(std::basic_stringbuf<charT,traits>::sync() == -1)
    return -1;

  basic_errorbuf<charT,traits>::str("");
  basic_errorbuf<charT,traits>::location("");
  return 0;
}

//=================================================================
//             Definition of basic error multi buffer
//-----------------------------------------------------------------
// Defines an error buffer class for pushing errors to an ostream
//=================================================================

template<class charT,class traits=std::char_traits<char> >
class basic_errormultibuf : public basic_errorbuf<charT,traits>
{
  public:
    typedef basic_errorbuf<charT,traits> base_type;

    typedef typename base_type::string_type         string_type;
    typedef typename base_type::error_priority_type error_priority_type;
    typedef typename base_type::error_info          error_info;

    explicit basic_errormultibuf();
    basic_errormultibuf(const basic_errormultibuf &sb);

    virtual ~basic_errormultibuf();

    void insert(const basic_errorstream<charT,traits> &o);
    void insert(std::ostream &o);
    void insert(std::ostream &o, const error_info &e);

    void restrict(error_restriction r, error_priority p);
    void restrict_range(error_priority p1, error_priority p2);

  protected:
    struct sink_info
    {
      typedef boost::shared_ptr<basic_errorbuf<charT,traits> > buf_ptr_type;
      error_restriction r;
      error_priority p1;
      error_priority p2;
      buf_ptr_type buf;

      explicit sink_info(buf_ptr_type b = buf_ptr_type())
      {
	p1 = no_priority;
	p2 = no_priority;
	r  = r_none;
	buf = b;
      }

      sink_info(const sink_info &s)
      {
        buf = buf_ptr_type();
        assign(s);
      }

      ~sink_info()
      {
      }

      sink_info &operator=(const sink_info &rhs)
      {
        assign(rhs);
	return *this;
      }

      void assign(const sink_info &rhs)
      {
	p1 = rhs.p1;
	p2 = rhs.p2;
	r  = rhs.r;

        buf = rhs.buf;
      }
    };

    typedef std::vector<sink_info> sink_vector;

    virtual int sync();
    virtual void push_error();
    sink_vector sinks;
};

//=================================================================
//     Implementation of basic error stream buffer (INLINE)
//=================================================================

template<class charT,class traits>
basic_errormultibuf<charT,traits>::basic_errormultibuf() :
        basic_errorbuf<charT,traits>()
{
  DEBUG(std::cout << "new basic_errormultibuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errormultibuf<charT,traits>::
basic_errormultibuf(const basic_errormultibuf<charT,traits> &sb) :
       basic_errorbuf<charT,traits>(sb), sinks(sb.sinks)
{
  DEBUG(std::cout << "new basic_errormultibuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errormultibuf<charT,traits>::~basic_errormultibuf()
{
  DEBUG(std::cout << "basic_errormultibuf being synced before delete at " << (void*) this << std::endl;)
  sync();
  DEBUG(std::cout << "deleting basic_errormultibuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::push_error()
{
  DEBUG(std::cout << "basic_errormultibuf::push_error() at " << (void*) this << std::endl;)
  if(!basic_errorbuf<charT,traits>::str().length() || !sinks.size()) return;

  for(typename sink_vector::iterator i = sinks.begin(); i != sinks.end(); ++i)
  {
    if(!i->buf.get()) continue;

    bool print = true;
    switch( i->r )
    {
      case r_eq: if( i->p1 != basic_errorbuf<charT,traits>::priority() ) print = false; break;
      case r_ne: if( i->p1 == basic_errorbuf<charT,traits>::priority() ) print = false; break;
      case r_lt: if( i->p1 <= basic_errorbuf<charT,traits>::priority() ) print = false; break;
      case r_gt: if( i->p1 >= basic_errorbuf<charT,traits>::priority() ) print = false; break;
      case r_le: if( i->p1 <  basic_errorbuf<charT,traits>::priority() ) print = false; break;
      case r_ge: if( i->p1 >  basic_errorbuf<charT,traits>::priority() ) print = false; break;

      case r_range:
        // range from [p1..p2] st p1<p2.  Don't print if p<p1 || p>p2
        if( i->p1 < i->p2 && (i->p1 > basic_errorbuf<charT,traits>::priority()
                               || i->p2 < basic_errorbuf<charT,traits>::priority()))
          print = false;
        else
        // range !(p2..p1) st p2<=p1.  Don't print if p2<p<p1
          if( i->p1 > basic_errorbuf<charT,traits>::priority() && i->p2 < basic_errorbuf<charT,traits>::priority())
            print = false;
        break;

      default: break;
    }

    if(!print)
      continue;

    i->buf->priority(basic_errorbuf<charT,traits>::priority());
    i->buf->filename(basic_errorbuf<charT,traits>::filename());
    i->buf->linenumber(basic_errorbuf<charT,traits>::linenumber());
    i->buf->str(basic_errorbuf<charT,traits>::str());
    if( basic_errorbuf<charT,traits>::raw_mode() )
      i->buf->set_raw_mode();
    else 
      i->buf->set_cooked_mode();
    i->buf->pubsync();
  }
}

template<class charT,class traits>
int basic_errormultibuf<charT,traits>::sync()
{
  DEBUG(std::cout << "basic_errormultibuf::sync() at " << (void*) this << std::endl;)
  if(!sinks.size())
    return -1;

  push_error();

  if(std::basic_stringbuf<charT,traits>::sync() == -1)
    return -1;

  basic_errorbuf<charT,traits>::str("");
  basic_errorbuf<charT,traits>::location("");
  return 0;
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::
  insert(const basic_errorstream<charT,traits> &o)
{
  DEBUG(std::cout << "basic_errormultibuf::insert1(" << (void*)o.rdbuf().get() << ") at " << (void*) this << std::endl;)
  sink_info s(o.rdbuf());
  sinks.push_back(s);
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::insert(std::ostream &o)
{
  DEBUG(std::cout << "basic_errormultibuf::insert2(ostream) at " << (void*) this << std::endl;)
  typename sink_info::buf_ptr_type t(new basic_errorstreambuf<charT,traits>(o));
  t->get_info() = basic_errorbuf<charT,traits>::get_info();
  sink_info s(t);
  sinks.push_back(s);
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::insert(std::ostream &o,
                                               const error_info &e)
{
  DEBUG(std::cout << "basic_errormultibuf::insert(ostream) at " << (void*) this << std::endl;)
  typename sink_info::buf_ptr_type t(new basic_errorstreambuf<charT,traits>(o));
  t->get_info() = e;
  sink_info s(t);
  sinks.push_back(s);
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::restrict(error_restriction r,
                                                 error_priority p)
{
  if(!sinks.size() || r == r_range) return;
  sinks.back().r = r;
  sinks.back().p1 = p;
  sinks.back().p2 = no_priority;
}

template<class charT,class traits>
void basic_errormultibuf<charT,traits>::restrict_range(error_priority p1,
                                                       error_priority p2)
{
  if(!sinks.size()) return;
  sinks.back().r = r_range;
  sinks.back().p1 = p1;
  sinks.back().p2 = p2;
}

}

#endif

