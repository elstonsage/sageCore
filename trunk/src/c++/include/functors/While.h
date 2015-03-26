#ifndef FUNCTORS_WHILE_H
#define FUNCTORS_WHILE_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "boost/bind.hpp"
#include "boost/mpl/not.hpp"
#include "boost/mpl/push_back.hpp"
#include "functors/ArgList.h"
#include "functors/make_function.h"
#include "functors/Misc.h"
#include "functors/const_ref.h"

namespace SAGE     {
namespace FUNCTORS { 

/// \brief Encapsulates a while(condtion) do { ... } control structure.
/// 
/// While is templatized on one argument:
///
/// ARGS - The input argument list (a boost::mpl::vector).
///
/// While stores two functors:
///
/// An continue condtion ("ExitFunctor") - Takes in ARGS and returns a bool indicating
/// whether or not to continue.
///
/// A body ("BodyFunctor") - Takes in ARGS and returns a bool indicating whether to
/// procede with the while structure or not (false is equivalent to a 'break' statement).
template<typename ARGS = boost::mpl::vector<> >
class While
{
public:

  /// @name Typedefs
  //@{

    DECLARE_ARGLIST;

    typedef void result_type;
    
    typedef typename make_function<bool, ARGS>::type ContinueFunctor;
    typedef typename make_function<bool, ARGS>::type BodyFunctor;
    
    typedef ReturnConstValue<bool, ARGS> DefaultContinueFunctor;
    typedef ReturnConstValue<bool, ARGS> DefaultBodyFunctor;

  //@}

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    While() : my_continue(DefaultContinueFunctor(true)), my_body(DefaultBodyFunctor(true)) { }
    
    ///
    /// Constructor #2.
    /// \param continue_fn The continue condition functor
    /// \param body_fn The body functor
    While(const ContinueFunctor & continue_fn, const BodyFunctor & body_fn) : my_continue(continue_fn), my_body(body_fn) { }
    
    ///
    /// Copy constructor.
    While(const While & other) : my_continue(other.my_continue), my_body(other.my_body) { }
    
    ///
    /// Operator=.
    While& operator=(const While & other) { if(this == &other) return *this; my_body = other.my_body; my_continue = other.my_continue; return *this; }
    
  //@}
  
  /// @name Setting the functors
  //@{
  
    ///
    /// Sets the loop body functor.
    void setBodyFunctor(const BodyFunctor & fn) { my_body = fn; }
  
    ///
    /// Sets the continue condition functor.
    void setContinueFunctor(const ContinueFunctor & fn) { my_continue = fn; }
  
  //@}
  
  /// @name Invoking the loop
  //@{
  
    void operator() () const
    {
      assertSeqSize<ARGS, 0> (); while(my_continue()) { if(!my_body()) break; }
    }
  
    void operator() (arg1_type a1) const
    {
      assertSeqSize<ARGS, 1> (); while(my_continue(a1)) { if(!my_body(a1)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<ARGS, 2> (); while(my_continue(a1, a2)) { if(!my_body(a1, a2)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<ARGS, 3> (); while(my_continue(a1, a2, a3)) { if(!my_body(a1, a2, a3)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<ARGS, 4> (); while(my_continue(a1, a2, a3, a4)) { if(!my_body(a1, a2, a3, a4)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<ARGS, 5> (); while(my_continue(a1, a2, a3, a4, a5)) { if(!my_body(a1, a2, a3, a4, a5)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<ARGS, 6> (); while(my_continue(a1, a2, a3, a4, a5, a6)) { if(!my_body(a1, a2, a3, a4, a5, a6)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<ARGS, 7> (); while(my_continue(a1, a2, a3, a4, a5, a6, a7)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<ARGS, 8> (); while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8)) break; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<ARGS, 9> (); while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8, a9)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8, a9)) break; }
    }

  //@}

private:

  ContinueFunctor my_continue;
  BodyFunctor     my_body;
};


} // End namespace FUNCTORS
} /// End namespace SAGE

#endif
