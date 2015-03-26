#ifndef FUNCTORS_MULTI_SELECTOR_H
#define FUNCTORS_MULTI_SELECTOR_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include "boost/bind.hpp"
#include "boost/mpl/not.hpp"
#include "boost/type_traits/remove_const.hpp"
#include "boost/type_traits/add_reference.hpp"
#include "functors/If.h"
#include "functors/ArgList.h"
#include "functors/const_ref.h"
#include "functors/is_void.h"

namespace SAGE {
namespace FUNCTORS {

/// Helper structs for the MultiSwitch
namespace MultiSwitch_Private {

///
/// Takes care of actually executing the selected functor and aggregating its result.
template<bool RETURN_IS_VOID> struct Executor { };

/// Specialized for non-void return type
template<> struct Executor<false>
{ 
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn());
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4));
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4, a5));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4, a5, a6));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4, a5, a6, a7));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7,
                                                                                                  typename MULTISWITCH::arg8_type a8)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4, a5, a6, a7, a8));
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7,
                                                                                                  typename MULTISWITCH::arg8_type a8,
                                                                                                  typename MULTISWITCH::arg9_type a9)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate(fn(a1, a2, a3, a4, a5, a6, a7, a8, a9));
  }
  
};

// Specialized for void:
template<> struct Executor<true>
{ 
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn)
  {
    fn();

    ((typename MULTISWITCH::Aggregator)ms).aggregate();
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1);
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                                      typename MULTISWITCH::arg2_type a2,
                                                                                                                      typename MULTISWITCH::arg3_type a3)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3);
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                                      typename MULTISWITCH::arg2_type a2,
                                                                                                                      typename MULTISWITCH::arg3_type a3,
                                                                                                                      typename MULTISWITCH::arg4_type a4)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4);
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                                      typename MULTISWITCH::arg2_type a2,
                                                                                                                      typename MULTISWITCH::arg3_type a3,
                                                                                                                      typename MULTISWITCH::arg4_type a4,
                                                                                                                      typename MULTISWITCH::arg5_type a5)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4, a5);
  }

  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4, a5, a6);
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4, a5, a6, a7);
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7,
                                                                                                  typename MULTISWITCH::arg8_type a8)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4, a5, a6, a7, a8);
  }
  
  template<typename MULTISWITCH, typename FN> static void go(MULTISWITCH & ms, const FN & fn, typename MULTISWITCH::arg1_type a1,
                                                                                                  typename MULTISWITCH::arg2_type a2,
                                                                                                  typename MULTISWITCH::arg3_type a3,
                                                                                                  typename MULTISWITCH::arg4_type a4,
                                                                                                  typename MULTISWITCH::arg5_type a5,
                                                                                                  typename MULTISWITCH::arg6_type a6,
                                                                                                  typename MULTISWITCH::arg7_type a7,
                                                                                                  typename MULTISWITCH::arg8_type a8,
                                                                                                  typename MULTISWITCH::arg9_type a9)
  {
    ((typename MULTISWITCH::Aggregator)ms).aggregate();
    
    fn(a1, a2, a3, a4, a5, a6, a7, a8, a9);
  }
  
};

//==========================
// DEFAULT AGGREGATOR STUFF
//==========================

template<typename T, bool RETURN_IS_VOID> class DefaultAggregatorBase;

template<typename T> class DefaultAggregatorBase<T, true>  { public: void aggregate ()                            const { } };
template<typename T> class DefaultAggregatorBase<T, false> { public: void aggregate (typename const_ref<T>::type) const { } };

template<typename I> class DefaultAggregator : public DefaultAggregatorBase<I, is_void<I>::value>
{
public:

  typedef void result_type;

  void reset() const { }

  result_type finalize() const { }
};

} // End namespace MultiSwitch_Private

/// \brief Given a input, selects zero or more conditional functors to operate on that input.
///
/// \par Introduction
///
/// Let's say you have an input object which can be processed according to a set of criteria. For each
/// criterion, if it is met (given the input), an action should be taken. If the criterion is not met,
/// no action is taken. Furthermore, the result of any actions taken has to be aggregated together somehow
/// in some kind of final report. This is what the MultiSwitch is for.
///
/// \par A quick example
///
/// We'll start with a simple example. The input in this case is an instance of the struct Person, which has
/// a int data member called "height", another int data member called "weight", and a std::string data member
/// called "name" (representing height, weight, and name, respectively).
///
/// For any Person instance, we want to do a few things. If the Person's height is over 50, we want to invoke
/// the function "wayTooTall". If the Person's weight is under 10, we want to invoke the function "tooLight".
/// And if the Person's name is "Lothario", we want to invoke "fatalErrorLotharioEncountered." We have three
/// sets of criteria, and for each set an action to take if the criteria is met.
///
/// Here's the code for this example. An explanation follows:
///
/// \code 
///
/// #include "functors/functors.h"
///
/// struct Person { int height; int weight; std::string name; }
///
/// bool isTooTall(const Person & person) { return person.height > 50; }
///
/// void wayTooTall(const Person & person) { std::cout << "way too Tall!!!" << std::endl; }
///
/// bool isTooLight(const Person & person) { return person.weight < 10; }
///
/// void tooLight(const Person & person) { std::cout << "you're too short." << std::endl; }
///
/// bool isLothario(const Person & person) { return person.name == "Lothario"; }
///
/// void fatalErrorLotharioEncountered(const Person & person) { exit(0); }
///
/// int main()
/// {
///   typedef SAGE::FUNCTORS::MultiSwitch<void, boost::mpl::vector<const Person &> > mtype;
///  
///   mtype m;
///  
///   m.addAction(boost::bind(&isTooTall,  _1), boost::bind(&wayTooTall,                    _1));
///   m.addAction(boost::bind(&isTooLight, _1), boost::bind(&tooLight,                      _1));
///   m.addAction(boost::bind(&isLothario, _1), boost::bind(&fatalErrorLotharioEncountered, _1));
///  
///   Person p;
///  
///   p.height = 55; p.weight = 45; p.name = "steve";
///  
///   m(p);
///  
///   p.height = 45; p.weight = 5; p.name = "Lothario";
///  
///   m(p);
///
///   return 0;
/// }
///
/// \endcode
///
/// First of all, let's take a look at the first line in the main() body:
///
/// \code typedef SAGE::FUNCTORS::MultiSwitch<void, boost::mpl::vector<const Person &> > mtype; \endcode
///
/// The first argument (void) specifies the return type of the MultiSwitch's "action"
/// functors. In this case, the MultiSwitch (and the action functors) will return void. The second argument
/// (boost::mpl::vector<const Person &>) is the input argument list. It indicates that the MultiSwitch, its
/// "if" functor, and its "action" functors, will all take a const reference to a Person as input.
///
/// Now, let's look at how the actions themselves are added. Adding an action to a MultiSwitch means specifying two things:
/// the conditions for executing the action, and the action to take if the conditions are met. All three of the
/// MultiSwitch::addAction() invocations give this information. The first line...
///
/// \code m.addAction(boost::bind(&isTooTall,  _1), boost::bind(&wayTooTall,                    _1)); \endcode
///
/// ...indicates the if isTooTall() returns true (given the input, of course), then the MultiSwitch should
/// execute wayTooTall().
///
/// After adding the actions, a Person instance is created. The person's attributes are assigned, and the MultiSwitch
/// is invoked (as a functor, that is). It in turn goes through its list of functor. For each pair of functors in the
/// list, it runs the first one to see if the given input meets the criteria for running the second. If it meets that
/// criteria (the "if" functor returns true) the MultiSwitch runs the "then" functor.
///
/// The first time the MultiSwitch is invoked in the above example, p's height is 55, p's weight is 45, and p's name is
/// "steve". The MultiSwitch will, then, invoke wayTooTall().
///
/// The second time the MultiSwitch is invoked, p's height is 45, p's weight is 5, and p's name is "Lothario". Accordingly,
/// the MultiSwitch will invoke tooLight() and fatalErrorLotharioEncountered().
/// 
/// \par Aggregators
///
/// There's one more (very important) aspect of a MultiSwitch to discuss: the Aggregator!
///
/// It may have occurred to you that it would be useful to process all the return values of the various "action"
/// functors. In our above example, the "action" functors have void return types. They could, however, have non-void
/// return types. Consider the following adaptation of the original algorithm:
///
/// If the Person's height is over 50, return security access "READ".
///
/// If the Person's weight is under 10, return security access "WRITE".
///
/// If the Person's name is "Lothario", return security access "EXECUTE".
///
/// Given a Person, figure out the Person's complete security access code (READ and/or WRITE and/or EXECUTE).
///
/// In the amended example, each "action" functor returns a security access code; these codes need to be "aggregated"
/// together somehow to represent a complete security access. The MultiSwitch is designed to let you accomplish this
/// kind of functionality.
///
/// Ok, so how is this done?
///
/// In addition to the return type and argument list, the MultiSwitch is actually templatized on a third argument:
/// AGGREGATOR. The AGGREGATOR, as a concept, has the following requirements:
///
/// 1. Default constructable
///
/// 2. Copy constructable
///
/// 3. Assignable
///
/// 4. Has a typedef called "result_type".
///
/// 5. Has a public function of the form: void reset() const. This function is invoked by the MultiSwitch
/// whenever its operator() is invoked.
///
/// 6. Has a public function of the form: void aggregate(const & R) const. This function is invoked by the MultiSwitch
/// whenever one of its "action" functors has been executed. The result of that execution is passed to the the
/// aggregate() function.
///
/// 7. Has a public function of the form: result_type finalize() const. This function is invoked by the MultiSwitch
/// whenever it finishes running its operator(). The MultiSwitch::operator() returns the value given by finalize().
///
/// 
template<typename R, typename ARGS, typename AGGREGATOR = MultiSwitch_Private::DefaultAggregator<R> >
class MultiSwitch : public AGGREGATOR
{

public:

  /// @name Typedefs
  //@{

    DECLARE_ARGLIST;

    typedef MultiSwitch_Private::Executor<is_void<R>::value> ExecutorType;
    typedef AGGREGATOR                                       BaseType;
    typedef AGGREGATOR                                       Aggregator;
    typedef typename Aggregator::result_type                 result_type;
    typedef If<R, args>                                      IfType;
    typedef std::vector<IfType>                              IfVector;

  //@}

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    MultiSwitch() { my_fns.clear(); }
    
    ///
    /// Copy constructor.
    MultiSwitch(const MultiSwitch & other) : BaseType(other), my_fns(other.my_fns) { }

    ///
    /// Assignment operator.
    MultiSwitch& operator=(const MultiSwitch & other) { if(this == &other) return *this; BaseType::operator=(other); my_fns = other.my_fns; return *this; }
    
  //@}

  /// @name Adding actions
  //@{

    ///
    /// Adds an "action" to this MultiSwitch.
    /// \param if_fn An "if" functor (takes ARGS as input and returns a boolean)
    /// \param action_fn An "action" functor (takes ARGS as input and returns R).
    void addAction(const typename IfType::PredicateFunctor & if_fn, const typename IfType::ActionFunctor & action_fn)
    {
      IfType c;
      c.setIfFunctor(if_fn);
      c.setThenFunctor(action_fn);

      my_fns.push_back(c);
    }
    
  //@}
  
  /// @name Invoking the selector
  //@{
  
    result_type operator() () const
    {
      assertSeqSize<args, 0> (); 
      
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4, a5); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4, a5, a6); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4, a5, a6, a7); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4, a5, a6, a7, a8); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      Aggregator::reset();

      for(typename IfVector::const_iterator s = my_fns.begin(); s != my_fns.end(); ++s)
        try { ExecutorType::go(*this, *s, a1, a2, a3, a4, a5, a6, a7, a8, a9); } catch(const std::exception &) { }
          
      return Aggregator::finalize();
    }
    
  //@}
  
private:

  IfVector my_fns;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
