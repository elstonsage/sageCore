#ifndef BUFFERED_ERROR_STREAM_H
#define BUFFERED_ERROR_STREAM_H

#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include "error/errorbuf.h"
#include "error/bufferederrorbuf.h"
#include "error/errorstream.h"

namespace SAGE {

template<class charT=char, class traits=std::char_traits<charT> >
class bufferederrorstream : public errorstream<charT, traits>
{
public:
    typedef errorstream<charT, traits> base_type;
    typedef charT char_type;
    typedef traits traits_type;
    typedef typename traits::pos_type pos_type;
    typedef typename traits::off_type off_type;
    typedef typename traits::int_type int_type;
    typedef std::char_traits<char_type> string_traits;
    typedef std::basic_string<char_type,string_traits,
                          std::allocator<char_type> > string_type;
    typedef basic_bufferederrorbuf<charT,traits> sb_type;
    typedef typename sb_type::error_info sb_error_info;
    typedef error_priority error_priority_type;

    typedef typename sb_type::iterator iterator;
    typedef typename sb_type::const_iterator const_iterator;
    typedef typename sb_type::error_type error_type;

    explicit bufferederrorstream(const basic_errorstream<charT,traits> &o);
    explicit bufferederrorstream(const bufferederrorstream<charT,traits> &s);

    inline iterator       begin();
    inline const_iterator begin() const;

    inline iterator       end();
    inline const_iterator end() const;

    inline error_type&       front();
    inline const error_type& front() const;

    inline error_type&       back();
    inline const error_type& back() const;

    inline void erase(iterator);

    inline void pop_front();
    inline void pop_back();

    inline void flush_buffer();
    inline void clear();
};

template<class charT,class traits>
bufferederrorstream<charT,traits>::
bufferederrorstream(const basic_errorstream<charT,traits> &o)
   : errorstream<charT,traits>(typename errorstream<charT,traits>::sb_type_ptr(new sb_type(o)))
{ }

template<class charT,class traits>
bufferederrorstream<charT,traits>::
bufferederrorstream(const bufferederrorstream<charT,traits> &o) :
          errorstream<charT,traits>(o.rdbuf()) { }

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::iterator
   bufferederrorstream<charT,traits>::begin()
{
  return static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->begin();
}

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::const_iterator
   bufferederrorstream<charT,traits>::begin() const
{
  return static_cast<const sb_type*>(errorstream<charT,traits>::sbuf.get())->begin();
}

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::iterator
   bufferederrorstream<charT,traits>::end()
{
  return static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->end();
}

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::const_iterator
   bufferederrorstream<charT,traits>::end() const
{
  return static_cast<const sb_type*>(errorstream<charT,traits>::sbuf.get())->end();
}

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::error_type&
   bufferederrorstream<charT,traits>::front()
{
  return static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->front();
}

template<class charT, class traits>
inline const typename bufferederrorstream<charT,traits>::error_type&
   bufferederrorstream<charT,traits>::front() const
{
  return static_cast<const sb_type*>(errorstream<charT,traits>::sbuf.get())->front();
}

template<class charT, class traits>
inline typename bufferederrorstream<charT,traits>::error_type&
   bufferederrorstream<charT,traits>::back()
{
  return static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->back();
}

template<class charT, class traits>
inline const typename bufferederrorstream<charT,traits>::error_type&
   bufferederrorstream<charT,traits>::back() const
{
  return static_cast<const sb_type*>(errorstream<charT,traits>::sbuf.get())->back();
}

template<class charT, class traits>
inline void
   bufferederrorstream<charT,traits>::erase(iterator i)
{
  static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->erase(i);
}

template<class charT, class traits>
inline void
   bufferederrorstream<charT,traits>::pop_front()
{
  static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->pop_front();
}

template<class charT, class traits>
inline void
   bufferederrorstream<charT,traits>::pop_back()
{
  static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->pop_back();
}

template<class charT, class traits>
inline void
   bufferederrorstream<charT,traits>::flush_buffer()
{
  static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->flush_buffer();
}

template<class charT, class traits>
inline void
   bufferederrorstream<charT,traits>::clear()
{
  static_cast<sb_type*>(errorstream<charT,traits>::sbuf.get())->clear();
}

}

#endif
