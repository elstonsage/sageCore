#ifndef FUNCTORS_FOR_H
#define FUNCTORS_FOR_H

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

/// \brief Encapsulates a loop structure as a standalone functor
///
/// This class lets you put together a conditional looping structure from your
/// own functors. A For is templatized on two arguments:
///
/// I - The typename of the index number. (For instance, in the loop for(int i = ...), 'i'
/// is the index number, and 'int' is the type of the index number. I is typedefed as index_type
/// in a For instance. The default value is 'size_t'.
///
/// ARGS - The list of input arguments to an invocation of For::operator(). This takes the
/// form of a boost::mpl::vector. The default value is an empty mpl vector.
///
/// There are four functors that For manages. By default, it uses ThrowException functors for all
/// four. Therefore, you will have to set all four functors in a For in order to be able to use it.
///
/// 1. Initialization ("InitFunctor"): Takes in ARGS and returns the initial value of the index variable.
///
/// 2. Continue condition ("ContinueFunctor"): Takes in ARGS + the index variable and returns a boolean indicating
/// whether or not the loop should continue.
///
/// 3. For body ("BodyFunctor"): Takes in ARGS + the index variable and returns a boolean (true means the
/// body exited normally; false indicates a 'break' statement to the loop).
///
/// 4. Index variable iteration ("IterFunctor"): Takes in ARGS + the index variable and returns the next
/// ("incremented") value of the index variable.
///
/// \par Example
///
/// Ok, so let's see a complete example. The following example shows how to create a For with an int as
/// the index value, going from values 0 to 10, and printing out the value for each iteration.
///
/// \code
///
/// int main()
/// {
///   SAGE::FUNCTORS::For<int> l;
///
///   l.setInitFunctor     (boost::lambda::constant(0));
///   l.setContinueFunctor (boost::lambda::_1 < 10));
///   l.setIterFunctor     (boost::lambda::_1 + 1);
///   l.setBodyFunctor     (std::cout << boost::lambda::_1);
///  
///   l();
///
///   return 0;
/// }
///
/// \endcode
///
template<typename I = size_t, typename ARGS = boost::mpl::vector<> >
class For
{
public:

  /// @name Typedefs
  //@{

    DECLARE_ARGLIST;

    typedef          void                                                                    result_type;
    typedef          I                                                                       index_type;
    typedef typename boost::mpl::push_back<ARGS, typename const_ref<index_type>::type>::type arglist_with_idx;
  
    typedef typename make_function<index_type, ARGS            >::type InitFunctor;
    typedef typename make_function<bool,       arglist_with_idx>::type BodyFunctor;
    typedef typename make_function<bool,       arglist_with_idx>::type ContinueFunctor;
    typedef typename make_function<index_type, arglist_with_idx>::type IterFunctor;
    
    typedef ThrowException<index_type, ARGS>             DefaultInitFunctor;
    typedef ThrowException<bool,       arglist_with_idx> DefaultBodyFunctor;
    typedef ThrowException<bool,       arglist_with_idx> DefaultContinueFunctor;
    typedef ThrowException<index_type, arglist_with_idx> DefaultIterFunctor;

  //@}

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    For() : my_init(DefaultInitFunctor()), my_body(DefaultBodyFunctor()), my_continue(DefaultContinueFunctor()), my_iter(DefaultIterFunctor()) { }
    
    ///
    /// Copy constructor.
    For(const For & other) : my_init(other.my_init), my_body(other.my_body), my_continue(other.my_continue), my_iter(other.my_iter) { }
    
    ///
    /// Operator=.
    For& operator=(const For & other) 
    { 
      if(this == &other) 
        return *this; 
      
      my_init = other.my_init; my_body = other.my_body; my_continue = other.my_continue; my_iter = other.my_iter; 
      
      return *this; 
    }
    
  //@}
  
  /// @name Setting the functors
  //@{
  
    ///
    /// Sets the initial condition functor.
    void setInitFunctor(const InitFunctor & fn) { my_init = fn; }

    ///
    /// Sets the loop body functor.
    void setBodyFunctor(const BodyFunctor & fn) { my_body = fn; }
  
    ///
    /// Sets the continue condition functor.
    void setContinueFunctor(const ContinueFunctor & fn) { my_continue = fn; }
  
    ///
    /// Sets the iteration functor.
    void setIterFunctor(const IterFunctor & fn) { my_iter = fn; }
  
  //@}
  
  /// @name Invoking the loop
  //@{
  
    void operator() () const
    {
      assertSeqSize<ARGS, 0> ();
      
      index_type i = my_init();
      
      while(my_continue(i))
      {
        my_body(i);
        i = my_iter(i);
      }
    }
  
    void operator() (arg1_type a1) const
    {
      assertSeqSize<ARGS, 1> ();
      
      index_type i = my_init(a1, i);
      
      while(my_continue(a1, i))
      {
        if(!my_body(a1, i)) break;
        i = my_iter(a1, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<ARGS, 2> ();
      
      index_type i = my_init(a1, a2, i);
      
      while(my_continue(a1, a2, i))
      {
        if(!my_body(a1, a2, i)) break;
        i = my_iter(a1, a2, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<ARGS, 3> ();
      
      index_type i = my_init(a1, a2, a3);
      
      while(my_continue(a1, a2, a3, i))
      {
        if(!my_body(a1, a2, a3, i)) break;
        i = my_iter(a1, a2, a3, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<ARGS, 4> ();
      
      index_type i = my_init(a1, a2, a3, a4);
      
      while(my_continue(a1, a2, a3, a4, i))
      {
        if(!my_body(a1, a2, a3, a4, i)) break;
        i = my_iter(a1, a2, a3, a4, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<ARGS, 5> ();
      
      index_type i = my_init(a1, a2, a3, a4, a5);
      
      while(my_continue(a1, a2, a3, a4, a5, i))
      {
        if(!my_body(a1, a2, a3, a4, a5, i)) break;
        i = my_iter(a1, a2, a3, a4, a5, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<ARGS, 6> ();
      
      index_type i = my_init(a1, a2, a3, a4, a5, a6);
      
      while(my_continue(a1, a2, a3, a4, a5, a6, i))
      {
        if(!my_body(a1, a2, a3, a4, a5, a6, i)) break;
        i = my_iter(a1, a2, a3, a4, a5, a6, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<ARGS, 7> ();
      
      index_type i = my_init(a1, a2, a3, a4, a5, a6, a7);
      
      while(my_continue(a1, a2, a3, a4, a5, a6, a7, i))
      {
        if(!my_body(a1, a2, a3, a4, a5, a6, a7, i)) break;
        i = my_iter(a1, a2, a3, a4, a5, a6, a7, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<ARGS, 8> ();
      
      index_type i = my_init(a1, a2, a3, a4, a5, a6, a7, a8);
      
      while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8, i))
      {
        if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8, i)) break;
        i = my_iter(a1, a2, a3, a4, a5, a6, a7, a8, i);
      }
    }
  
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<ARGS, 9> ();
      
      index_type i = my_init(a1, a2, a3, a4, a5, a6, a7, a8, a9);
      
      while(my_continue(a1, a2, a3, a4, a5, a6, a7, a8, a9, i))
      {
        if(!my_body(a1, a2, a3, a4, a5, a6, a7, a8, a9, i)) break;
        i = my_iter(a1, a2, a3, a4, a5, a6, a7, a8, a9, i);
      }
    }

  //@}

private:

  InitFunctor     my_init;
  BodyFunctor     my_body;
  ContinueFunctor my_continue;
  IterFunctor     my_iter;
};


} // End namespace FUNCTORS
} /// End namespace SAGE

#endif
