#ifndef ERROR_STREAM_H
#define ERROR_STREAM_H

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <iostream>
#include "error/errorbuf.h"

namespace SAGE {

struct ErrorNoInit { };

template<class charT=char, class traits=std::char_traits<charT> >
class basic_errorstream : public std::basic_ostream<charT, traits>
{
public:
    typedef std::basic_ostream<charT,traits>	base_type;

    typedef typename base_type::char_type    char_type;
    typedef typename base_type::traits_type  traits_type;
    typedef typename base_type::pos_type     pos_type;
    typedef typename base_type::off_type     off_type;
    typedef typename base_type::int_type     int_type;

    typedef std::char_traits<char_type> string_traits;
    typedef std::basic_string<char_type,string_traits> string_type;
    typedef basic_errorbuf<charT,traits> sb_type;
    typedef boost::shared_ptr<sb_type>   sb_type_ptr;
    typedef typename sb_type::error_info sb_error_info;
    typedef error_priority error_priority_type;

    explicit basic_errorstream(const sb_type_ptr& sb = sb_type_ptr());
    explicit basic_errorstream(const basic_errorstream &s);
    virtual ~basic_errorstream();
    
    basic_errorstream &operator=(const basic_errorstream &s);

    //lint -e{1511} We want to override the basic_ostream rdbuf
    inline sb_type_ptr rdbuf() const;

    //lint -e{1511} We want to override the basic_ostream rdbuf
    void rdbuf(const sb_type_ptr& s);

    inline sb_error_info &get_info();
    inline string_type str() const;
    inline void str(const string_type& str);
    inline void prefix(const string_type &pre);
    inline string_type prefix() const;
    inline void suffix(const string_type &pre);
    inline string_type suffix() const;
    inline void location(const string_type &file, size_t line = ((size_t)-1));
    inline void filename(const string_type &f);
    inline string_type filename() const;
    inline void linenumber(size_t p);
    inline size_t linenumber() const;
    inline void priority(error_priority_type p);
    inline error_priority_type priority() const;
    inline void offset(size_t p);
    inline size_t offset() const;
    inline void indent(size_t p);
    inline size_t indent() const;
    inline void line_width(size_t p);
    inline size_t line_width() const;
    inline bool   raw_mode() const;
    inline void   set_raw_mode();
    inline bool   cooked_mode() const;
    inline void   set_cooked_mode();

protected:
    sb_type_ptr sbuf;
};

template<class charT,class traits>
inline
basic_errorstream<charT,traits>::basic_errorstream(const sb_type_ptr& sbuffer) :
      std::basic_ostream<charT,traits>(sbuffer.get()), sbuf(sbuffer)
{
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)sbuffer.get() << ")"<< std::endl;)
  rdbuf(rdbuf());
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)rdbuf().get() << ")"<< std::endl;)
}

//lint -e{1738} We don't use the copy constructor of basic_ostream because, at
//              least under KAI C++, there are instantiation problems with basic_ios.
template<class charT,class traits>
basic_errorstream<charT,traits>::
basic_errorstream(const basic_errorstream<charT,traits> &s) :
      std::basic_ostream<charT,traits>(s.sbuf.get()), sbuf(s.sbuf)
{
  DEBUG(std::cout << "basic_errorstream constructor(sb=" << (void*)sb << ")"<< std::endl;)
}

template<class charT,class traits>
basic_errorstream<charT,traits>::~basic_errorstream()
{
  DEBUG(std::cout << "deleting basic_errorstream at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_errorstream<charT,traits> &
basic_errorstream<charT,traits>::operator=
     (const basic_errorstream<charT,traits> &s)
{
  if(&s != this)
    rdbuf(s.rdbuf());

  return *this;
}

template<class charT,class traits>
void basic_errorstream<charT,traits>::rdbuf(const sb_type_ptr& sbuffer)
{
  if(rdbuf() == sbuffer && sbuf == sbuffer) return;

  DEBUG(std::cout << "basic_errorstream::rdbuf(sb=" << (void*)sbuffer.get() << ")"<< std::endl;)
  sbuf = sbuffer;

  //lint -e{534} Ok to ignore return value.
  std::basic_ostream<charT,traits>::rdbuf(sbuffer.get());
}

template<class charT,class traits>
inline
typename basic_errorstream<charT,traits>::sb_type_ptr 
basic_errorstream<charT,traits>::rdbuf() const 
{ return sbuf; }

//lint -e{36} Incomprehensible, and likely spurious, error
template<class charT, class traits, class T>
basic_errorstream<char,traits>& operator<<(basic_errorstream<char,traits>& out, const T& t)
{
  static_cast<std::basic_ostream<charT,traits>&>(out) << t;
  return out;
}

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline
typename basic_errorstream<charT,traits>::sb_error_info & 
basic_errorstream<charT,traits>::get_info() 
{ return rdbuf()->get_info(); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type 
basic_errorstream<charT,traits>::str() const 
{ return rdbuf()->str(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::str(const string_type& str_arg)
{ rdbuf()->str(str_arg); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline 
void basic_errorstream<charT,traits>::prefix(const string_type& str_arg)
{ rdbuf()->prefix(str_arg); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type 
basic_errorstream<charT,traits>::prefix() const
{ return rdbuf()->prefix(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline 
void basic_errorstream<charT,traits>::suffix(const string_type& str_arg)
{ rdbuf()->suffix(str_arg); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type
basic_errorstream<charT,traits>::suffix() const
{ return rdbuf()->suffix(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline 
void basic_errorstream<charT,traits>::location(const string_type& file, 
                                                     size_t line)
{ rdbuf()->location(file,line); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline void basic_errorstream<charT,traits>::filename(const string_type &f)
{ rdbuf()->filename(f); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::string_type
basic_errorstream<charT,traits>::filename() const
{return rdbuf()->filename(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::linenumber(size_t l)
{ rdbuf()->linenumber(l); }

template<class charT,class traits>
inline size_t basic_errorstream<charT,traits>::linenumber() const
{ return rdbuf()->linenumber(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline void basic_errorstream<charT,traits>::priority(error_priority_type pr)
{ rdbuf()->priority(pr); }

template<class charT,class traits> 
inline typename basic_errorstream<charT,traits>::error_priority_type
basic_errorstream<charT,traits>::priority() const
{ return rdbuf()->priority(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits> 
inline void basic_errorstream<charT,traits>::offset(size_t o)
{ rdbuf()->offset(o); }

template<class charT,class traits>
inline size_t basic_errorstream<charT,traits>::offset() const
{ return rdbuf()->offset(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::indent(size_t i)
{ rdbuf()->indent(i); }

template<class charT,class traits>
inline size_t basic_errorstream<charT,traits>::indent() const
{ return rdbuf()->indent(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::line_width(size_t i)
{ rdbuf()->line_width(i); }

template<class charT,class traits>
inline size_t basic_errorstream<charT,traits>::line_width() const
{ return rdbuf()->line_width(); }

template<class charT,class traits>
inline bool basic_errorstream<charT,traits>::raw_mode() const
{ return rdbuf()->raw_mode(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::set_raw_mode()
{ rdbuf()->set_raw_mode(); }

template<class charT,class traits>
inline bool basic_errorstream<charT,traits>::cooked_mode() const
{ return rdbuf()->cooked_mode(); }

//lint -e{1762} Function isn't const since it changes internal state
template<class charT,class traits>
inline void basic_errorstream<charT,traits>::set_cooked_mode()
{ rdbuf()->set_cooked_mode(); }

template<class charT=char, class traits=std::char_traits<charT> >
class errorstream : public basic_errorstream<charT, traits>
{
public:
    //lint -e{1516} Override expected
    typedef basic_errorstream<charT,traits> base_type;

    typedef typename base_type::char_type           char_type;
    typedef typename base_type::traits_type         traits_type;
    typedef typename base_type::pos_type            pos_type;
    typedef typename base_type::off_type            off_type;
    typedef typename base_type::int_type            int_type;

    typedef typename base_type::string_traits       string_traits;
    typedef typename base_type::string_type         string_type;
    typedef typename base_type::error_priority_type error_priority_type;

    typedef basic_errorstreambuf<charT,traits> errstream_type;
    typedef basic_errorbuf<charT,traits>       base_sb_type;
    typedef boost::shared_ptr<base_sb_type>    base_sb_type_ptr;

    explicit errorstream(const ErrorNoInit& );
    explicit errorstream(std::ostream &out = std::cerr);
    explicit errorstream(std::ostream &out, const string_type& str);
    errorstream(const errorstream<charT,traits> &s);


    virtual ~errorstream();

protected:
    errorstream(const base_sb_type_ptr& sb);
};

template<class charT,class traits>
errorstream<charT,traits>::errorstream(const ErrorNoInit& ) : 
              basic_errorstream<charT,traits>( basic_errorstream<charT,traits>::rdbuf() ) { }

//lint -e{1732} new is ok without assignment since we're using shared_ptr
template<class charT,class traits>
errorstream<charT,traits>::errorstream(std::ostream &ostr) : 
    basic_errorstream<charT,traits>
        (typename basic_errorstream<charT,traits>::sb_type_ptr
            (new errstream_type(ostr)))
{ }

template<class charT,class traits>
errorstream<charT,traits>::
errorstream(const base_sb_type_ptr& sbuffer) : 
              basic_errorstream<charT,traits>(sbuffer) 
{ 
  DEBUG((std::cout << "errorstream non-init constructor(sb=" << (void*)basic_errorbuf<charT,traits>::rdbuf().get() << ")"<< std::endl;))
}

template<class charT,class traits>
errorstream<charT,traits>::errorstream
   (std::ostream&      ostr,
    const string_type& st)
   : basic_errorstream<charT,traits>(sb_type_ptr(new errstream_type(ostr,st)))
{ }

template<class charT,class traits>
errorstream<charT,traits>::
errorstream(const errorstream<charT,traits> &o) :
      basic_errorstream<charT,traits>(o) { }

template<class charT,class traits>
errorstream<charT,traits>::~errorstream()
{ }

template<class charT=char, class traits=std::char_traits<charT> >
class errormultistream : public errorstream<charT, traits>
{
public:
    //lint -e{1516} Override expected
    typedef errorstream<charT,traits> base_type;

    typedef typename base_type::char_type            char_type;
    typedef typename base_type::traits_type          traits_type;
    typedef typename base_type::pos_type             pos_type;
    typedef typename base_type::off_type             off_type;
    typedef typename base_type::int_type             int_type;

    typedef typename base_type::string_traits        string_traits;
    typedef typename base_type::string_type          string_type;
    typedef typename base_type::error_priority_type  error_priority_type;
    typedef typename base_type::sb_error_info        sb_error_info;
    
    typedef basic_errormultibuf<charT,traits> multibuf_type;
    typedef boost::shared_ptr<multibuf_type> multibuf_type_ptr;

    errormultistream();
    errormultistream(const errormultistream<charT,traits> &s);
    virtual ~errormultistream();

    inline void insert(const basic_errorstream<charT,traits> &o);
    inline void insert(std::ostream &o);
    inline void insert(std::ostream &o, const sb_error_info &e);
    inline void restrict(error_restriction r, error_priority p);
    inline void restrict_range(error_priority p1, error_priority p2);
};

template<class charT,class traits>
errormultistream<charT,traits>::
errormultistream() : 
          errorstream<charT,traits>(typename errorstream<charT,traits>::base_sb_type_ptr(new multibuf_type())) { }

template<class charT,class traits>
errormultistream<charT,traits>::
errormultistream(const errormultistream<charT,traits> &o) :
          errorstream<charT,traits>(o.rdbuf()) { }

template<class charT,class traits>
errormultistream<charT,traits>::
~errormultistream()
{ }

template<class charT,class traits> 
inline void 
errormultistream<charT,traits>::
  insert(const basic_errorstream<charT,traits> &o)
{ static_cast<multibuf_type*>(basic_errorstream<charT,traits>::rdbuf().get())->insert(o); }

template<class charT,class traits>
inline 
void
errormultistream<charT,traits>::
  insert(std::ostream &o)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->insert(o); }

template<class charT,class traits>
inline 
void
errormultistream<charT,traits>::
  insert(std::ostream &o, const sb_error_info &e)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->insert(o,e); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::
  restrict(error_restriction r, error_priority p)
{ static_cast<multibuf_type*>(basic_errorstream<charT,traits>::rdbuf().get())->restrict(r,p); }

template<class charT,class traits>
inline void
errormultistream<charT,traits>::
  restrict_range(error_priority p1, error_priority p2)
{ static_cast<multibuf_type*>(basic_errorbuf<charT,traits>::rdbuf().get())->restrict_range(p1,p2); }

typedef errorstream<char> cerrorstream; 
typedef errormultistream<char> cerrormultistream; 

extern cerrorstream sage_cerr;
extern cerrorstream sage_clog;
extern cerrorstream sage_cout;

class errorstream_init 
{
  static long count;
public:
  errorstream_init();
 ~errorstream_init();
};

//lint -e{1502} Warning about no non-static members supressed
//static errorstream_init _do_errorstream_init;

void manual_errorstream_init();
  
}

#endif

