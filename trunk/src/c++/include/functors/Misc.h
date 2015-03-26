#ifndef FUNCTORS_MISC_H
#define FUNCTORS_MISC_H

#include "boost/mpl/vector.hpp"
#include "functors/ArgList.h"

namespace SAGE     {
namespace FUNCTORS {

/// \brief Functor that will throw an exception whenever it's invoked
///
template<typename R, typename ARGS, typename E = std::exception> class ThrowException
{
  public:
  
  DECLARE_ARGLIST;

  R operator() ()                                                                                                  const { throw E(); } 
  R operator() (arg1_type)                                                                                         const { throw E(); } 
  R operator() (arg1_type, arg2_type)                                                                              const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type)                                                                   const { throw E(); }   
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type)                                                        const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type)                                             const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type)                                  const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type)                       const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type)            const { throw E(); } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type, arg9_type) const { throw E(); } 
};

/// \brief Functor that will return a single const value whenever it is invoked
/// 
template<typename R, typename ARGS = boost::mpl::vector<> > class ReturnConstValue
{
public:

  DECLARE_ARGLIST;

  ReturnConstValue() { }
  
  ReturnConstValue(const R & r) : my_val(r) { }

  ReturnConstValue(const ReturnConstValue & other) : my_val(other.my_val) { }

  ReturnConstValue& operator=(const ReturnConstValue & other) { if(this == &other) return *this; my_val = other.my_val; return *this; }

  R operator() ()                                                                                                  const { return my_val; } 
  R operator() (arg1_type)                                                                                         const { return my_val; } 
  R operator() (arg1_type, arg2_type)                                                                              const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type)                                                                   const { return my_val; }   
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type)                                                        const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type)                                             const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type)                                  const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type)                       const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type)            const { return my_val; } 
  R operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type, arg9_type) const { return my_val; } 

private:

  R my_val;
};

/// \brief Functor with void return type that will do nothing whenever it is invoked
///
template<typename ARGS = boost::mpl::vector<> > class DoNothing
{
public:

  DECLARE_ARGLIST;

  void operator() ()                                                                                                  const { } 
  void operator() (arg1_type)                                                                                         const { } 
  void operator() (arg1_type, arg2_type)                                                                              const { } 
  void operator() (arg1_type, arg2_type, arg3_type)                                                                   const { }   
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type)                                                        const { } 
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type)                                             const { } 
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type)                                  const { } 
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type)                       const { } 
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type)            const { } 
  void operator() (arg1_type, arg2_type, arg3_type, arg4_type, arg5_type, arg6_type, arg7_type, arg8_type, arg9_type) const { } 
};


} // End namespace FUNCTORS
} // End namespace SAGE

#endif
