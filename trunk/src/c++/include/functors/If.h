#ifndef FUNCTORS_CONDITIONAL_H
#define FUNCTORS_CONDITIONAL_H

#include "boost/mpl/not.hpp"
#include "functors/is_void.h"
#include "functors/const_ref.h"
#include "functors/make_not_void.h"
#include "functors/make_function.h"
#include "functors/Misc.h"

namespace SAGE     {
namespace FUNCTORS {

/// \brief Encapsulates a conditional control structure (if-then-else)
///
/// \par Introduction
/// 
/// Let's say you've got an if-then-else sequence (such as "if name is 'steve',
/// exit program; otherwise, print the name to the screen). Each step in the
/// conditional sequence (if, then, and else) can be expressed as a standalone
/// functor. For instance:
///
/// For the given input "name":
///
/// 'If' functor: Return true if name 'steve'.
///
/// 'Then' functor: Exit program.
///
/// 'Else' functor: Print name.
///
/// The 'If' class lets you put together the functors above and use the
/// If instance as a functor on its own.
///
/// \par Template arguments
///
/// A If instance is templatized on return type and an argument list. The return
/// type is what the If::operator() will return; the argument list (a boost::mpl vector)
/// lists the argument types that will be passed to the If::operator().
///
/// \par The 'if' functor
///
/// To set the 'if' functor, use the setIfFunctor() function. The functor must take as input the arglist
/// given by ARGS (see above), and return a boolean.
///
/// \par The 'then' and 'else' functors
///
/// The 'then' and 'else' functors must take ARGS and return R (see above).
///
/// \par Default behavior
///
/// By default, a If instance sets its three functors (if, then, and else)
/// to ThrowException functors. That is, if you create a If instance but assign
/// it no functors, it will throw exceptions whenever used.
///
/// \par A complete example
///
/// Let's say you've got the following if-then-else algorithm: Given age, if age > 20, exit the program; 
/// other, print the age to the screen. The three distinct steps -- if age > 20, exiting the program, and
/// printing the age to the screen -- must exist somewhere in your program as "functorizable" components.
/// In the following example, I've shown how to implement the conditional structure entirely with non-member
/// functions:
///
/// \code
///
/// bool ageOver20(int age) { return age  > 20; }
///
/// void exitNow(int age) { exit(0); }
///
/// void reportAge(int age) { std::cout << age; }
///
/// int main()
/// {
///   typedef SAGE::FUNCTORS::If<void, boost::mpl::vector<int> > ctype;
///  
///   ctype c;
///  
///   c.setIfFunctor   (boost::bind(&ageOver20, _1));
///   c.setThenFunctor (boost::bind(&exitNow,   _1));
///   c.setElseFunctor (boost::bind(&reportAge, _1));
///  
///   c(15);
///   c(25);
///
///   return 0;
/// }
///
/// \endcode
///
template<typename R = void, typename ARGS = boost::mpl::vector<> >
class If
{
public:

  /// @name Typedefs
  //@{

    DECLARE_ARGLIST;

    typedef R                                        result_type;
    typedef typename make_function<bool, ARGS>::type PredicateFunctor;
    typedef typename make_function<R, ARGS>::type    ActionFunctor;

  //@}

  /// @name Constructors
  //@{

    ///
    /// Default constructor. Initializes all functors to ThrowException.
    If() : my_if(ThrowException<bool, args> ()), my_then(ThrowException<result_type, args> ()), my_else(ThrowException<result_type, args> ()) { }

    ///
    /// Constructor
    /// \param if_fn The 'if' functor
    /// \param then_fn The 'then' functor
    /// \param else_fn The 'else' functor
    If(const PredicateFunctor & if_fn, const ActionFunctor & then_fn, const ActionFunctor & else_fn) : my_if(if_fn), my_then(then_fn), my_else(else_fn) { }

    ///
    /// Constructor
    /// \param if_fn The 'if' functor
    /// \param r1 A constant value to return if the condition if the 'if' functor is met.
    /// \param r2 A constant value to return if the condition if the 'if' functor is not met.
    If(
      const PredicateFunctor & if_fn, 
      typename const_ref<typename make_not_void<R>::type>::type r1, 
      typename const_ref<typename make_not_void<R>::type>::type r2) 
      : 
      my_if   (if_fn), 
      my_then (ReturnConstValue<typename make_not_void<R>::type, ARGS>(r1)), 
      my_else (ReturnConstValue<typename make_not_void<R>::type, ARGS>(r2)) 
      
    { 
      // Makes sure that the user hasn't tried to use this constructor if the return type is void.
      BOOST_STATIC_ASSERT(( boost::mpl::not_<is_void<R> >::value ));
    }
    
    ///
    /// Copy constructor.
    If(const If & other) : my_if(other.my_if), my_then(other.my_then), my_else(other.my_else) { }
    
    ///
    /// Assignment operator.
    If& operator=(const If & other) { if(this == &other) return *this; my_if = other.my_if; my_then = other.my_then; my_else = other.my_else; return *this; }
    
  //@}
  
  /// @name Setting the functors
  //@{
  
    ///
    /// Assigns the 'if' functor (which takes ARGS and returns a boolean).
    void setIfFunctor(const PredicateFunctor & if_fn) { my_if = if_fn; }

    ///
    /// Assigns the 'then' functor (action to take if 'if' criteria is met).
    void setThenFunctor(const ActionFunctor & then_fn) { my_then = then_fn; }

    ///
    /// Assigns the 'else' functor (action to take if 'if' criteria is NOT met).
    void setElseFunctor(const ActionFunctor & else_fn) { my_else = else_fn; }
  
  //@}
  
  /// @name Invoking the functor
  //@{
  
    result_type operator() () const
    {
      assertSeqSize<args, 0> ();
      
      return my_if() ? my_then() : my_else();
    }
  
    result_type operator() (arg1_type a1) const
    {
      assertSeqSize<args, 1> ();
      
      return my_if(a1) ? my_then(a1) : my_else(a1);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<args, 2> ();
      
      return my_if(a1, a2) ? my_then(a1, a2) : my_else(a1, a2);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<args, 3> ();
      
      return my_if(a1, a2, a3) ? my_then(a1, a2, a3) : my_else(a1, a2, a3);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<args, 4> ();
      
      return my_if(a1, a2, a3, a4) ? my_then(a1, a2, a3, a4) : my_else(a1, a2, a3, a4);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<args, 5> ();
      
      return my_if(a1, a2, a3, a4, a5) ? my_then(a1, a2, a3, a4, a5) : my_else(a1, a2, a3, a4, a5);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<args, 6> ();
      
      return my_if(a1, a2, a3, a4, a5, a6) ? my_then(a1, a2, a3, a4, a5, a6) : my_else(a1, a2, a3, a4, a5, a6);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<args, 7> ();
      
      return my_if(a1, a2, a3, a4, a5, a6, a7) ? my_then(a1, a2, a3, a4, a5, a6, a7) : my_else(a1, a2, a3, a4, a5, a6, a7);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<args, 8> ();
      
      return my_if(a1, a2, a3, a4, a5, a6, a7, a8) ? my_then(a1, a2, a3, a4, a5, a6, a7, a8) : my_else(a1, a2, a3, a4, a5, a6, a7, a8);
    }
  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<args, 9> ();
      
      return my_if(a1, a2, a3, a4, a5, a6, a7, a8, a9) ? my_then(a1, a2, a3, a4, a5, a6, a7, a8, a9) : my_else(a1, a2, a3, a4, a5, a6, a7, a8, a9);
    }
  
  //@}

private:

  PredicateFunctor my_if;
  ActionFunctor    my_then;
  ActionFunctor    my_else;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
