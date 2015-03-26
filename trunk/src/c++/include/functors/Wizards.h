#ifndef FUNCTORS_WIZARDS_H
#define FUNCTORS_WIZARDS_H

#include "functors/While.h"
#include "functors/For.h"
#include "functors/If.h"
#include "boost/lambda/lambda.hpp"

namespace SAGE     {
namespace FUNCTORS {

///
/// Returns an If<...> instance via C++-like syntax.
///
/// Syntax:
///
/// if_(condition) [ then_action ] .else_ [ else_action ]
///
/// Example:
///
/// \code
///
/// SAGE::FUNCTORS::If<> i = 
///   SAGE::FUNCTORS::if_<> (boost::lambda::constant(true))
///     [ std::cout << boost::lambda::constant(1) ]
///     .else_
///     [ std::cout << boost::lambda::constant(2) ];
///
/// \endcode
///
template<typename R = void, typename ARGS = boost::mpl::vector<> >
struct if_
{
  typedef If<R, ARGS> IfType;

  struct IfTypeHolder
  {
    IfType operator[] (const typename IfType::ActionFunctor & else_fn)
    {
      i.setElseFunctor(else_fn);
      return i;
    }

    IfType i;
  };

  explicit if_(const typename IfType::PredicateFunctor & predicate)
  {
    else_.i.setIfFunctor(predicate);
  }
  
  if_ operator[] (const typename IfType::ActionFunctor & then_fn)
  {
    else_.i.setThenFunctor(then_fn);
    return *this;
  }
  
  operator IfType () const { return else_.i; }

  IfTypeHolder else_;
};

///
/// Returns a While<...> instance via C++-like syntax.
///
/// Syntax:
///
/// while_(continue condition) [ body ];
///
/// For example:
///
/// \code
/// While<> w = while_<> (ReturnConstValue<bool>(true)) [ ReturnConstValue<bool>(false) ];
/// w();
/// \endcode
///
template<typename ARGS = boost::mpl::vector<> >
struct while_
{
  typedef While<ARGS> WhileType;
  
  explicit while_(const typename WhileType::ContinueFunctor & continue_fn)
  {
    w.setContinueFunctor(continue_fn);
  }

  WhileType operator[] (const typename WhileType::BodyFunctor & body)  
  {
    w.setBodyFunctor(body);
    return w;
  }

  WhileType w;
};

///
/// Returns a DoWhile<...> instance via C++-like syntax.
///
/// Syntax:
///
/// do_ <...> () [ body ] .while_ ( continue condition );
///
/// For example:
///
/// \code
///
/// DoWhile<> w = do_<> () [ std::cout << boost::lambda::constant('i') ] .while (boost::lambda::constant(true));
///
/// \endcode
///
template<typename ARGS = boost::mpl::vector<> >
struct do_
{
  typedef DoWhile<ARGS> DoWhileType;
  
  struct DoWhileHolder
  {
    DoWhileType while_(const typename DoWhileType::ContinueFunctor & continue_fn)
    {
      d.setContinueFunctor(continue_fn);
      return d;
    }
    
    DoWhileType d;
  };
  
  DoWhileHolder operator[] (const typename DoWhileType::BodyFunctor & body)
  {
    DoWhileHolder hldr;
    hldr.d.setBodyFunctor(body);
    return hldr;
  }
};

///
/// Returns a For<...> instance via C++-like syntax.
///
/// Syntax:
///
/// for_(initial condition, continue condition, iteration functor) [ body ]
///
/// Example:
///
/// \code
///
/// int main()
/// {
///   SAGE::FUNCTORS::For<> l = 
///     SAGE::FUNCTORS::for_<>(boost::lambda::constant(0), boost::lambda::_1 >= 10, boost::lambda::_1 + 1)
///     [
///       std::cout << boost::lambda::_1
///     ];
///  
///   l();
///
///   return 0;
/// }
///
/// \endcode
template<typename I = size_t, typename ARGS = boost::mpl::vector<> >
struct for_
{
  typedef For<I, ARGS> ForType;
  
  for_(const typename ForType::InitFunctor     & init_,
       const typename ForType::ContinueFunctor & continue_,
       const typename ForType::IterFunctor     & iter_)
  {
    f.setInitFunctor     (init_);
    f.setContinueFunctor (continue_);
    f.setIterFunctor     (iter_);
  }
  
  ForType operator[] (const typename ForType::BodyFunctor & body)
  {
    f.setBodyFunctor(body);
    return f;
  }
  
  ForType f;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
