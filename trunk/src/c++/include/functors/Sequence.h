#ifndef FUNCTORS_SEQUENCE_H
#define FUNCTORS_SEQUENCE_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "boost/bind.hpp"
#include "boost/mpl/not.hpp"
#include "functors/make_function.h"
#include "functors/Misc.h"
#include "functors/const_ref.h"

namespace SAGE     {
namespace FUNCTORS { 

/// \brief Manages an ordered sequence of functors
///
/// \par Prerequisites
///
/// First of all, make sure your familiar with Functor's before you get started
/// learning this class. Consider it a prerequisite.
///
/// \par Template arguments reference
///
/// R - The return type of the functors that make up the sequence
///
/// ARGS - A boost::mpl::vector<> specifying which arguments the sequence functors take
///
/// \par Introduction
///
/// Let's say you've got some procedure that can be expressed as a sequence of
/// logical steps. For instance, consider the calculation of BMI,
/// (where BMI = (weight / (height^2)) * 703 ). We could decompose this equation
/// into discrete steps like this:
///
/// 1. BMI = weight
///
/// 2. BMI /= (height^2)
///
/// 3. BMI *= 703
///
/// The Sequence class exists to organize logical steps (coded as functors) to
/// carry out some procedure. The basic idea is that you program each step as a 
/// standalone functor, then add each functor in the correct procedural order
/// to your Sequence. Your Sequence can then be called simply (with operator() )
/// to calculate the final result.
///
/// Let's consider another example before going on to the details of how to use
/// a Sequence. The previous example was mathematical in nature; what about a logical
/// sequence that is non-mathematical? Consider the following procedure for
/// an assembly-line building of a car:
///
/// 1. Requisition chassis
///
/// 2. Install engine
///
/// 3. Install transmission
///
/// 4. Add body
///
/// 5. Test for acceptable performance
///
/// 6. Deliver to buyer
///
/// Each one of the above steps is a discrete, logical step in the overall process of
/// assembling a car. The Sequence class is ideally suited to help you organize those
/// steps.
///
/// \par A quick example to get you started
///
/// When you create a Sequence, it is templatized on the same types on which its
/// constituent Functor's are templatized. That is, a Sequence conatins a list of Functors
/// of the \b same type. The Sequence itself is designed to work as functor as well; it
/// has the complete standard Functor interface / typedefs / etc.
///
/// When the Sequence is invoked (via operator() ), the Sequence will in turn invoke
/// each of its constituent Functors sequentially. It will pass to those Functor's the
/// input argument(s) you have passed to the Sequence's operator(). It will return the
/// returned value of the \b last functor in the list. (Actually, it's a bit more complicated
/// than that, but we'll save that for later).
///
/// Ok, so let's see how to code the earlier BMI example from above. Recall that there
/// are three discrete steps (set BMI equal to weight; divide it by height squared;
/// multiply it by 703). The first thing we'll need is three functors, each of which
/// has a pointer to the BMI variable (note that the functors store their results
/// by modifiying pointers, \b not by returning calculated values. More on this later).
///
/// \code
/// 
/// double BMI = 0.0;
/// 
/// class SetToWeight 
/// { 
/// public:
///   double * bmi;
///   explicit SetToWeight(double * _bmi) : bmi(_bmi) { }
///   void operator() (double weight, double height) const { *bmi = weight; }
/// };
/// 
/// class DivideByHeightSquared
/// {
/// public:
///   double * bmi;
///   explicit DivideByHeightSquared(double * _bmi) : bmi(_bmi) { }
///   void operator() (double weight, double height) const { *bmi /= (height * height); }
/// };
/// 
/// class MultiplyBy703
/// {
/// public:
///   double * bmi;
///   explicit MultiplyBy703(double * _bmi) : bmi(_bmi) { }
///   void operator() (double weight, double height) const { *bmi *= 703; }
/// };
///
/// \endcode
/// 
/// Now we just need to add a main function in which we create a Sequence, add the functors
/// we want, and invoke the Sequence a few times for practice.
///
/// \code
///
/// #include "functors/Sequence.h"
/// 
/// int main()
/// {
///   SAGE::Sequence<void, boost::mpl::vector<double, double> > seq;
///   
///   seq.addFunctor("set to weight",      SetToWeight           (&BMI));
///   seq.addFunctor("divide by height^2", DivideByHeightSquared (&BMI));
///   seq.addFunctor("multiply by 703",    MultiplyBy703         (&BMI));
///   
///   seq(5, 10);
///   
///   std::cout << BMI; // BMI has been modified by the previous statement; let's see what it is...
///   
///   seq(15, 100);
///   
///   std::cout << BMI; // BMI has been modified by the previous statement; let's see what it is...
///   
///   return 0;
/// }
///
/// \endcode
///
/// What happens when seq(x, x) is invoked? The Sequence invokes, one by one, each of its functors
/// with the given input. It then returns the return value of the last step. In the above example,
/// our functors' return type is void, so the Sequence doesn't return anything.
///
/// \par Step-wise evaluation (what if a step fails?!)
///
/// What if a step has "failed" and the sequence should abort? Does the Sequence allow for
/// this circumstance?
///
/// Yes, it does! When you instantiate your Sequence, you can optionally pass an Evaluator
/// object. This object must obey the following functor requirements: it must take as input a const
/// reference to each Functor's output (the return type of the Sequence), and it must return a boolean
/// indicating whether or not this output was valid. For instance, if each step in your Sequence returns
/// an int, and an int less than 10 should cause the sequence to abort, the following evaluator functor
/// would work:
///
/// \code
/// class LessThan10 { public: bool operator() (int i) const { return i < 10; } };
/// \endcode
///
/// Unless you specify your own functor when you instantiate your Sequence, your Sequence will
/// assume every step is always successful.
///
///
template<typename R, typename ARGS = boost::mpl::vector<> >
class Sequence
{
public:

  // Typedefs for the template parameters
  
  typedef R result_type;

  // Typedef for the stored functor type:

  typedef typename make_function<R, ARGS>::type FunctorType;
  
  // Typedefs for arguments:

  DECLARE_ARGLIST;

  // Evaluator stuff:

  typedef ReturnConstValue<bool, boost::mpl::vector<typename const_ref<result_type>::type> > DefaultEvaluatorType;

  typedef typename make_function<bool, boost::mpl::vector<typename const_ref<result_type>::type> >::type EvaluatorFunctor;

  // Miscellaneous:

  /// \brief Represents a single "step" in the sequence (a name, evaluator, and action)
  struct Step
  {
    Step() { }
    
    Step(const Step & other) : name(other.name), evaluator(other.evaluator), functor(other.functor) { }

    Step& operator=(const Step & other) { if(this == &other) return *this; name = other.name; evaluator = other.evaluator; functor = other.functor; return *this; }
    
    std::string      name;
    EvaluatorFunctor evaluator;
    FunctorType      functor;
  };

  typedef std::vector<Step> StepVector;

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    Sequence()
    {
      my_steps   . clear();
      my_success = false;
      my_fn_idx  = 0;
    }
    
    ///
    /// Copy constructor.
    Sequence(const Sequence & other) : 
      my_steps    (other.my_steps), 
      my_success  (other.my_success),
      my_fn_idx   (other.my_fn_idx)
    { }

    ///
    /// Assignment operator.
    Sequence& operator=(const Sequence & other)
    {
      if(this == &other)
        return *this;
        
      my_steps   = other.my_steps;
      my_success = other.my_success;
      my_fn_idx  = other.my_fn_idx;
      
      return *this;
    }
    
  //@}

  /// @name Adding / querying functors
  //@{

    ///
    /// Adds another step to the sequence.
    /// NOTE: This version of addStep pairs the given 'step' with a default evaluator (which always
    /// returns true).
    /// \param fn The step to add (must be a functor taking ARGS and returning R)
    /// \param name The name to assign to the functor.
    void addStep(const FunctorType & fn, const std::string & name)
    {
      Step step;
      
      step.name      = name;
      step.evaluator = DefaultEvaluatorType(true);
      step.functor   = fn;
      
      my_steps.push_back(step);
    }

    ///
    /// Adds another step to the sequence.
    /// \param fn The step to add (must be a functor taking ARGS and returning R)
    /// \param eval The evaluator to add (must be a functor taking a const R & and returning a boolean).
    /// \param name The name to assign to the functor.
    void addStep(const FunctorType & fn, const EvaluatorFunctor & eval, const std::string & name)
    {
      Step step;
      
      step.name      = name;
      step.evaluator = eval;
      step.functor   = fn;
      
      my_steps.push_back(step);
    }

  //@}
  
  /// @name Sequence results
  //@{
  
    ///
    /// Returns whether or not the most recent operator() invocation was successful
    /// (according to the evaluator).
    bool getSuccess() const { return my_success; }
    
    ///
    /// Returns the index number of the last executed functor (from the most recent
    /// invocation of operator() ).
    const Step & getLastExecutedStep() const { return my_steps[my_fn_idx]; }
    
  //@}

  /// @name Invoking the sequence
  //@{

    result_type operator() () const
    {
      assertSeqSize<args, 0> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor());
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor());
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  
    
    result_type operator() (arg1_type a1) const
    {
      assertSeqSize<args, 1> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<args, 2> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<args, 3> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<args, 4> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<args, 5> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<args, 6> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<args, 7> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<args, 8> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<args, 9> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size() - 1; ++my_fn_idx)
      {
        result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8, a9));
        
        my_success = my_steps[my_fn_idx].evaluator(r);
       
        if(!getSuccess())
          return r;
      }
      
      result_type r(my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8, a9));
        
      my_success = my_steps[my_fn_idx].evaluator(r);
       
      return r;
    }  

  //@}
  
private:

          StepVector my_steps;
  mutable bool       my_success;
  mutable size_t     my_fn_idx;
};

/// Specialized for a void return type
template<typename ARGS>
class Sequence<void, ARGS>
{
public:

  // Typedefs for the template parameters
  
  typedef void result_type;

  // Typedef for the stored functor type:

  typedef typename make_function<result_type, ARGS>::type FunctorType;
  
  // Typedefs for arguments:

  DECLARE_ARGLIST;

  // Evaluator stuff:

  typedef ReturnConstValue<bool, boost::mpl::vector< > > DefaultEvaluatorType;

  typedef typename make_function<bool, boost::mpl::vector< > >::type EvaluatorFunctor;

  // Miscellaneous:

  struct Step
  {
    Step() { }
    Step(const Step & other) : name(other.name), evaluator(other.evaluator), functor(other.functor) { }
    Step& operator=(const Step & other) { if(this == &other) return *this; name = other.name; evaluator = other.evaluator; functor = other.functor; return *this; }
    
    std::string      name;
    EvaluatorFunctor evaluator;
    FunctorType      functor;
  };

  typedef std::vector<Step> StepVector;

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    Sequence()
    {
      my_steps   . clear();
      my_success = false;
      my_fn_idx  = 0;
    }
    
    ///
    /// Copy constructor.
    Sequence(const Sequence & other) : 
      my_steps    (other.my_steps), 
      my_success  (other.my_success),
      my_fn_idx   (other.my_fn_idx)
    { }

    ///
    /// Assignment operator.
    Sequence& operator=(const Sequence & other)
    {
      my_steps   = other.my_steps;
      my_success = other.my_success;
      my_fn_idx  = other.my_fn_idx;
      
      return *this;
    }
    
  //@}

  /// @name Adding / querying functors
  //@{

    ///
    /// Adds another step to the sequence.
    /// NOTE: This version of addStep pairs the given 'step' with a default evaluator (which always
    /// returns true).
    /// \param fn The step to add (must be a functor taking ARGS and returning R)
    /// \param name The name to assign to the functor.
    template<typename FN>
    void addStep(const FN & fn, const std::string & name)
    {
      Step step;
      
      step.name      = name;
      step.evaluator = DefaultEvaluatorType(true);
      step.functor   = fn;
      
      my_steps.push_back(step);
    }

    ///
    /// Adds another step to the sequence.
    /// \param fn The step to add (must be a functor taking ARGS and returning R)
    /// \param eval The evaluator to add (must be a functor taking a const R & and returning a boolean).
    /// \param name The name to assign to the functor.
    template<typename FN, typename EVAL>
    void addStep(const FN & fn, const EVAL & eval, const std::string & name)
    {
      Step step;
      
      step.name      = name;
      step.evaluator = eval;
      step.functor   = fn;
      
      my_steps.push_back(step);
    }

  //@}
  
  /// @name Sequence results
  //@{
  
    ///
    /// Returns whether or not the most recent operator() invocation was successful
    /// (according to the evaluator).
    bool getSuccess() const { return my_success; }
    
    ///
    /// Returns the index number of the last executed functor (from the most recent
    /// invocation of operator() ).
    const Step & getLastExecutedStep() const { return my_steps[my_fn_idx]; }
    
  //@}

  /// @name Invoking the sequence
  //@{

    result_type operator() () const
    {
      assertSeqSize<args, 0> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor();
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  
    
    result_type operator() (arg1_type a1) const
    {
      assertSeqSize<args, 1> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<args, 2> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<args, 3> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<args, 4> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  
    
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<args, 5> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<args, 6> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<args, 7> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<args, 8> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  

    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<args, 9> ();
      assert(my_steps.size());

      for(my_fn_idx = 0; my_fn_idx < my_steps.size(); ++my_fn_idx)
      {
        my_steps[my_fn_idx].functor(a1, a2, a3, a4, a5, a6, a7, a8, a9);
        
        my_success = my_steps[my_fn_idx].evaluator();
       
        if(!getSuccess())
          return;
      }
    }  

  //@}
  
private:

          StepVector         my_steps;
  mutable bool               my_success;
  mutable size_t             my_fn_idx;
};


} // End namespace FUNCTORS
} /// End namespace SAGE

#endif
