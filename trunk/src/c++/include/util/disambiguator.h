#ifndef DISAMBIGUATOR_H
#define DISAMBIGUATOR_H

#include "boost/type_traits/add_reference.hpp"
#include "boost/type_traits/add_const.hpp"

#define DISAMBIGUATE(CLASSNAME, DATATYPE) \
class CLASSNAME                                               \
{                                                             \
  public:                                                     \
                                                              \
    typedef boost::add_reference<DATATYPE>::type Reference;   \
                                                              \
    typedef boost::add_const<Reference>::type ConstReference; \
                                                              \
    typedef CLASSNAME self;                                   \
                                                              \
    explicit CLASSNAME(Reference d)                           \
        : my_data(d) { }                                      \
                                                              \
    CLASSNAME(const self& o)                                  \
        : my_data(o.my_data) { }                              \
                                                              \
    operator Reference () { return my_data; }                 \
                                                              \
    operator ConstReference () const { return my_data; }      \
                                                              \
    Reference operator() () { return my_data; }               \
                                                              \
    ConstReference operator() () const { return my_data; }    \
                                                              \
  private:                                                    \
                                                              \
    Reference my_data;                                        \
};

#endif

