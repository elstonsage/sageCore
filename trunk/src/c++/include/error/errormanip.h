#ifndef ERROR_MANIP_H
#define ERROR_MANIP_H

//=====================================================================
//  File:       errormanip.h
//
//  Purpose:    This header file defines error stream io manipulators
//              for streams derived from basic_errorstream.
//
//  Author:     Kevin Jacobs
//
//  History:    Version 0.10
//
//  Copyright (c) 1998 R.C. Elston.  All Rights Reserved
//=====================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include "error/errorstream.h"

namespace SAGE {

#define error_location location(__FILE__,__LINE__)

#if 0
template<class charT, class traits>
basic_errorstream<charT, traits>& raw_mode(basic_errorstream<charT, traits>& es)
{
  es.raw_mode();
  return es;
}

template<class charT, class traits>
basic_errorstream<charT, traits>& cooked_mode(basic_errorstream<charT, traits>& es)
{
  es.cooked_mode();
  return es;
}
#endif

template <class T, class charT, class traits = std::char_traits<charT> >
class basic_erroromanip 
{

#if WHY_IS_THIS_BROKEN
    friend
    basic_errorstream<charT, traits>&
    operator<< (basic_errorstream<charT, traits>& os, 
                const basic_erroromanip<T, charT, traits>& a);
#endif

public :
 
    typedef  charT               char_type;
    typedef  traits              traits_type;

    typedef  typename traits::pos_type    pos_type;
    typedef  typename traits::off_type    off_type;
    typedef  typename traits::int_type    int_type;
 
protected:
 
    typedef  basic_errorstream<charT, traits>    ios_type;
    typedef  ios_type&  (* pf_type) (ios_type&, T);
 
public :

    inline basic_erroromanip (pf_type, T);

#ifdef WHY_IS_THIS_PROKEN
private:
#endif
    pf_type    pf;
    T          manarg;
};

template <class T, class charT, class traits>
inline
basic_erroromanip<T, charT, traits>::
basic_erroromanip (pf_type pf_arg, T manarg_arg)
: pf (pf_arg), manarg (manarg_arg) { }       

template <class T>
class error_omanip : public basic_erroromanip<T, char, std::char_traits<char> >
{
private:

    typedef  basic_errorstream<char, std::char_traits<char> > ios_type;
    typedef  ios_type&  (* pf_type) (ios_type&, T);

public :

    inline error_omanip (pf_type pf_arg, T arg)
       : basic_erroromanip<T, char, std::char_traits<char> > (pf_arg, arg)
   { }

};                      

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
fpriority (basic_errorstream<charT, traits>& stream, error_priority p)
{
    // set priority
    stream.priority(p);
    return stream;
}

inline
error_omanip<error_priority>
priority(error_priority p)
{
    return error_omanip<error_priority>(&fpriority, p);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
fprefix (basic_errorstream<charT, traits>& stream, const std::string &p)
{
    // set prefix
    stream.prefix(p);
    return stream;
}

inline
error_omanip<const std::string &>
prefix(const std::string &p)
{
    return error_omanip<const std::string &>(&fprefix, p);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
fsuffix (basic_errorstream<charT, traits>& stream, 
         const std::string &p)
{
    // set suffix
    stream.suffix(p);
    return stream;
}

inline
error_omanip<const std::string &>
suffix(const std::string &p)
{
    return error_omanip<const std::string &>(&fsuffix, p);
}

struct error_loc
{
  error_loc(const std::string f, size_t l) : file(f), line(l) {}
  const std::string file;
  size_t line;
};

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
flocation (basic_errorstream<charT, traits>& stream, error_loc loc)
{
    // set location
    stream.location(loc.file, loc.line);
    return stream;
}

inline
error_omanip<error_loc>
location(const std::string &file, size_t line = ((size_t)-1))
{
    return error_omanip<error_loc>(&flocation, error_loc(file,line));
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
ffilename (basic_errorstream<charT, traits>& stream, 
         const std::string &p)
{
    // set filename
    stream.filename(p);
    return stream;
}

inline
error_omanip<const std::string &>
filename(const std::string &p)
{
    return error_omanip<const std::string &>(&ffilename, p);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
flinenumber (basic_errorstream<charT, traits>& stream,  size_t i)
{
    // set linenumber
    stream.linenumber(i);
    return stream;
}

inline
error_omanip<size_t>
linenumber(size_t i)
{
    return error_omanip<size_t>(&flinenumber, i);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
findent (basic_errorstream<charT, traits>& stream,  size_t i)
{
    // set indent
    stream.indent(i);
    return stream;
}

inline
error_omanip<size_t>
indent(size_t i)
{
    return error_omanip<size_t>(&findent, i);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
foffset (basic_errorstream<charT, traits>& stream,  size_t o)
{
    // set offset
    stream.offset(o);
    return stream;
}

inline
error_omanip<size_t>
offset(size_t o)
{
    return error_omanip<size_t>(&foffset, o);
}

template <class charT, class traits>
inline
basic_errorstream<charT, traits>&
fline_width (basic_errorstream<charT, traits>& stream,  size_t o)
{
    // set line_width
    stream.line_width(o);
    return stream;
}

inline
error_omanip<size_t>
line_width(size_t o)
{
    return error_omanip<size_t>(&fline_width, o);
}

template <class T, class charT, class traits>
inline
basic_errorstream<charT, traits>&
operator<< (basic_errorstream<charT, traits>& os,
            const basic_erroromanip<T, charT, traits>& a)
{
    (*a.pf) (os, a.manarg);
    return os;
}                

}

#endif
