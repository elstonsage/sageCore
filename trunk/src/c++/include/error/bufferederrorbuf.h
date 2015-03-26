#ifndef BUFFERED_ERROR_BUF_H
#define BUFFERED_ERROR_BUF_H

//==================================================================
//  File:       bufferederrorbuf.h
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
#include <list>
#include "error/errorbuf.h"

#ifndef DEBUG
#define DEBUG(x)
#endif

namespace SAGE {

//=================================================================
//             Definition of basic buffered error buffer             
//-----------------------------------------------------------------
// Defines an error buffer class for pushing errors to an ostream 
//=================================================================

template<class charT,class traits=std::char_traits<char> >
class basic_bufferederrorbuf : public basic_errorbuf<charT,traits> 
{
  public:
    typedef basic_errorbuf<charT,traits> base_type;
    typedef std::basic_string<charT,traits> string_type;
    typedef basic_errorbuf<charT,traits> err_buf;
    typedef typename basic_errorbuf<charT,traits>::error_info error_info;
    typedef error_priority error_priority_type;
    typedef basic_bufferederrorbuf<charT,traits> buf_type;

    struct error_type
    {
      error_priority      priority;
      string_type         filename;
      size_t              linenumber;
      string_type         str;
    };

    typedef std::list<error_type>               error_list;
    typedef typename error_list::iterator       iterator;
    typedef typename error_list::const_iterator const_iterator;

    explicit basic_bufferederrorbuf(const basic_errorstream<charT,traits>& sb);
    explicit basic_bufferederrorbuf(const basic_errorstream<charT,traits>& sb, bool flush);
    explicit basic_bufferederrorbuf(const basic_bufferederrorbuf &sb);
    explicit basic_bufferederrorbuf(const basic_bufferederrorbuf &sb, bool flush);
    virtual ~basic_bufferederrorbuf();

// Add list operators

    iterator       begin();
    const_iterator begin() const;

    iterator       end();
    const_iterator end() const;

    error_type&       front();
    const error_type& front() const;

    error_type&       back();
    const error_type& back() const;

    void erase(iterator);

    void pop_front();
    void pop_back();

    void flush_buffer();
    void clear();

  protected:

    virtual int sync();
    virtual void push_error();

    error_list my_errors;

    boost::shared_ptr<basic_errorbuf<charT, traits> > my_buf;

    bool flush_on_destroy;
};

//=================================================================
//     Implementation of basic buffered error buffer (INLINE)       
//=================================================================

template<class charT,class traits>
basic_bufferederrorbuf<charT,traits>::
    basic_bufferederrorbuf(const basic_errorstream<charT,traits>& sb) :
        basic_errorbuf<charT,traits>(), flush_on_destroy(true)
{
  my_buf = sb.rdbuf();
  DEBUG(std::cout << "new basic_bufferederrorbuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_bufferederrorbuf<charT,traits>::
    basic_bufferederrorbuf(const basic_bufferederrorbuf<charT,traits>& sb) :
        basic_errorbuf<charT,traits>(), flush_on_destroy(sb.flush_on_destroy)
{
  my_buf = sb;
  DEBUG(std::cout << "new basic_bufferederrorbuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_bufferederrorbuf<charT,traits>::
    basic_bufferederrorbuf(const basic_errorstream<charT,traits>& sb, bool f) :
        basic_errorbuf<charT,traits>(), flush_on_destroy(f)
{
  my_buf = sb.rdbuf();

  DEBUG(std::cout << "new basic_bufferederrorbuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_bufferederrorbuf<charT,traits>::
basic_bufferederrorbuf(const basic_bufferederrorbuf<charT,traits> &sb, bool f) :
       basic_errorbuf<charT,traits>(sb), flush_on_destroy(f)
{ 
  my_buf = sb.rdbuf();

  DEBUG(std::cout << "new basic_bufferederrorbuf at " << (void*) this << std::endl;)
}

template<class charT,class traits>
basic_bufferederrorbuf<charT,traits>::~basic_bufferederrorbuf() 
{ 
  if(flush_on_destroy) sync(); 

  DEBUG(std::cout << "deleting basic_bufferederrorbuf at " << (void*) this << std::endl;)
}

template<class charT,class traits> 
void basic_bufferederrorbuf<charT,traits>::push_error()
{
  DEBUG(std::cout << "basic_bufferederrorbuf::push_error() at " << (void*) this << std::endl;)

  if(!basic_errorbuf<charT,traits>::str().length()) return;

  my_errors.push_back(error_type());
  
  error_type& e = back();
  
  e.priority   = basic_errorbuf<charT,traits>::priority();
  e.filename   = basic_errorbuf<charT,traits>::filename();
  e.linenumber = basic_errorbuf<charT,traits>::linenumber();
  e.str        = basic_errorbuf<charT,traits>::str();
}

template<class charT,class traits> 
int basic_bufferederrorbuf<charT,traits>::sync() 
{
  DEBUG(std::cout << "basic_bufferederrorbuf::sync() at " << (void*) this << std::endl;)

  push_error();

  if(std::basic_stringbuf<charT,traits>::sync() == -1)
    return -1;

  basic_errorbuf<charT,traits>::str("");
  basic_errorbuf<charT,traits>::location("");
  return 0;
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::iterator
    basic_bufferederrorbuf<charT,traits>::begin()
{
  return my_errors.begin();
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::const_iterator
    basic_bufferederrorbuf<charT,traits>::begin() const
{
  return my_errors.begin();
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::iterator
    basic_bufferederrorbuf<charT,traits>::end()
{
  return my_errors.end();
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::const_iterator
    basic_bufferederrorbuf<charT,traits>::end() const
{
  return my_errors.end();
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::error_type&
    basic_bufferederrorbuf<charT,traits>::front()
{
  return my_errors.front();
}

template<class charT,class traits> 
inline const typename basic_bufferederrorbuf<charT,traits>::error_type&
    basic_bufferederrorbuf<charT,traits>::front() const
{
  return my_errors.front();
}

template<class charT,class traits> 
inline typename basic_bufferederrorbuf<charT,traits>::error_type&
    basic_bufferederrorbuf<charT,traits>::back()
{
  return my_errors.back();
}

template<class charT,class traits> 
inline const typename basic_bufferederrorbuf<charT,traits>::error_type&
    basic_bufferederrorbuf<charT,traits>::back() const
{
  return my_errors.back();
}


template<class charT,class traits> 
inline void basic_bufferederrorbuf<charT,traits>::erase(iterator i)
{
  my_errors.erase(i);
}

template<class charT,class traits> 
inline void basic_bufferederrorbuf<charT,traits>::pop_front()
{
    my_errors.pop_front();
}

template<class charT,class traits> 
inline void basic_bufferederrorbuf<charT,traits>::pop_back()
{
    my_errors.pop_back();
}

template<class charT,class traits> 
inline void basic_bufferederrorbuf<charT,traits>::flush_buffer()
{
  while(my_errors.size())
  {
    my_buf->priority(front().priority);
    my_buf->filename(front().filename);
    my_buf->linenumber(front().linenumber);
    my_buf->str(front().str);

    my_buf->pubsync();
    
    pop_front();
  }

  clear();
}

template<class charT,class traits> 
inline void basic_bufferederrorbuf<charT,traits>::clear()
{
  my_errors.clear();
}

}

#endif
