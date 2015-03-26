#ifndef FUNCTORS_DO_WHILE_H
#define FUNCTORS_DO_WHILE_H

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

template<typename ARGS = boost::mpl::vector<> >
class DoWhile
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
    DoWhile() : my_continue(DefaultContinueFunctor(true)), my_body(DefaultBodyFunctor(true)) { }
    
    ///
    /// Constructor #2.
    /// \param continue_fn The continue condition functor
    /// \param body_fn The body functor
    DoWhile(const ContinueFunctor & continue_fn, const BodyFunctor & body_fn) : my_continue(continue_fn), my_body(body_fn) { }
    
    ///
    /// Copy constructor.
    DoWhile(const DoWhile & other) : my_continue(other.my_continue), my_body(other.my_body) { }
    
    ///
    /// Operator=.
    DoWhile& operator=(const DoWhile & other) { if(this == &other) return *this; my_body = other.my_body; my_continue = other.my_continue; return *this; }
    
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
      assertSeqSize<ARGS, 0> (); if(!my_body()) return; while(my_continue()) { if(!my_body()) return; }
    }
  
    void operator() (arg1_type a1) const
    {
      assertSeqSize<ARGS, 1> (); if(!my_body(a1)) return; while(my_continue(a1)) { if(!my_body(a1)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<ARGS, 2> (); if(!my_body(a1, a2)) return; while(my_continue(a1, a2)) { if(!my_body(a1, a2)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<ARGS, 3> (); if(!my_body(a1, a2, a3)) return; while(my_continue(a1, a2, a3)) { if(!my_body(a1, a2, a3)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<ARGS, 4> (); if(!my_body(a1, a2, a3, a4)) return; while(my_continue(a1, a2, a3, a4)) { if(!my_body(a1, a2, a3, a4)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<ARGS, 5> (); if(!my_body(a1, a2, a3, a4, a5)) return; while(my_continue(a1, a2, a3, a4, a5)) { if(!my_body(a1, a2, a3, a4, a5)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<ARGS, 6> (); if(!my_body(a1, a2, a3, a4, a5, a6)) return; while(my_continue(a1, a2, a3, a4, a5, a6)) { if(!my_body(a1, a2, a3, a4, a5, a6)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<ARGS, 7> (); if(!my_body(a1, a2, a3, a4, a5, a6, a7)) return; while(my_continue(a1, a2, a3, a4, a5, a6, a7)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<ARGS, 8> (); if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8)) return; while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8)) return; }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<ARGS, 9> (); if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8, a9)) return; while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8, a9)) { if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8, a9)) return; }
    }

  //@}

private:

  ContinueFunctor my_continue;
  BodyFunctor my_body;
};


} /// End namespace FUNCTORS
} /// End namespace SAGE

#endif
