#ifndef FUNCTORS_SWITCH_H
#define FUNCTORS_SWITCH_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include "boost/bind.hpp"
#include "boost/mpl/not.hpp"
#include "boost/type_traits/remove_const.hpp"
#include "boost/type_traits/add_reference.hpp"
#include "functors/Misc.h"
#include "functors/make_function.h"

namespace SAGE     {
namespace FUNCTORS {

/// \brief Given a classifiable input, chooses a preset processor to handle that input.
///
/// \par Introduction
///
/// Let's say you've got an input object that can be classified according to some classification
/// scheme. Furthermore, you want to choose an action on the basis of the that classification.
/// That's what the Switch class is for!
///
/// \par A quick example
///
/// 
template<typename C, typename R, typename ARGS = boost::mpl::vector<> >
class Switch
{

public:

  /// @name Typedefs
  //@{

    DECLARE_ARGLIST;

    typedef R                                                  result_type;
    typedef C                                                  Classification;
    typedef typename make_function<Classification, args>::type ClassifierFunctor;
    typedef typename make_function<result_type, args>::type    FunctorType;
    typedef std::map<Classification, FunctorType>              FunctorMap;

  //@}

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    Switch() : my_default_fn(ThrowException<R, ARGS> ()), my_classifier(ThrowException<Classification, ARGS> ()) { my_fns.clear(); }
    
    ///
    /// Copy constructor.
    Switch(const Switch & other) : my_default_fn(other.my_default_fn), my_fns(other.my_fns), my_classifier(other.my_classifier) { }

    ///
    /// Assignment operator.
    Switch& operator=(const Switch & other) { if(this == &other) return *this; my_fns = other.my_fns; my_default_fn = other.my_default_fn; my_classifier = other.my_classifier; return *this; }
    
  //@}

  /// @name Adding / removing functors
  //@{

    ///
    /// Sets the given functor to process the given classification.
    /// \param c The classification id that this functor will handle
    /// \param f The functor to process the input args
    void setFunctor(Classification c, const FunctorType & f) { my_fns[c] = FunctorType(f); }
    
    ///
    /// Sets the functor for handling any input set that (1) unclassifiable, or (2) lacks a
    /// classification-specific functor.
    void setDefaultFunctor(const FunctorType & f) { my_default_fn = FunctorType(f); }
    
    ///
    /// Sets the classifier functor.
    void setClassifier(const ClassifierFunctor & f) { my_classifier = f;  }
    
  //@}
  
  /// @name Accessing functors
  //@{

    ///
    /// Returns a const reference to the vector of functors.
    const FunctorMap & getFunctors() const { return my_fns; }
    
    ///
    /// Returns a const reference to the vector of functors.
    FunctorMap & getFunctors() { return my_fns; }
    
    ///
    /// Locates the functor for the given classification index.
    /// If no functor for that index was entered, it returns the default functor.
    const FunctorType & locateFunctor(Classification c) const
    {
      typename FunctorMap::const_iterator fn = my_fns.find(c);
      
      if(fn == my_fns.end()) return my_default_fn;
      else                   return fn->second;
    }
    
    ///
    /// Returns a const reference to the default functor.
    const FunctorType & getDefaultFunctor() const { return my_default_fn; }
    
    ///
    /// Returns a const reference to the classifier functor.
    const ClassifierFunctor & getClassifier() const { return my_classifier; }

  //@}
  
  /// @name Invoking the selector
  //@{
  
    void operator() () const
    {
      assertSeqSize<args, 0> (); try { locateFunctor(my_classifier())(); } catch(const std::exception & e) { my_default_fn(); }
    }
    
    void operator() (arg1_type a1) const
    {
      assertSeqSize<args, 1> (); try { locateFunctor(my_classifier(a1))(a1); } catch(const std::exception & e) { my_default_fn(a1); }
    }
    
    void operator() (arg1_type a1, arg2_type a2) const
    {
      assertSeqSize<args, 2> (); try { locateFunctor(my_classifier(a1, a2))(a1, a2); } catch(const std::exception & e) { my_default_fn(a1, a2); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3) const
    {
      assertSeqSize<args, 3> (); try { locateFunctor(my_classifier(a1, a2, a3))(a1, a2, a3); } catch(const std::exception & e) { my_default_fn(a1, a2, a3); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4) const
    {
      assertSeqSize<args, 4> (); try { locateFunctor(my_classifier(a1, a2, a3, a4))(a1, a2, a3, a4); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
    {
      assertSeqSize<args, 5> (); try { locateFunctor(my_classifier(a1, a2, a3, a4, a5))(a1, a2, a3, a4, a5); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4, a5); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
    {
      assertSeqSize<args, 6> (); try { locateFunctor(my_classifier(a1, a2, a3, a4, a5, a6))(a1, a2, a3, a4, a5, a6); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4, a5, a6); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
    {
      assertSeqSize<args, 7> (); try { locateFunctor(my_classifier(a1, a2, a3, a4, a5, a6, a7))(a1, a2, a3, a4, a5, a6, a7); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4, a5, a6, a7); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
    {
      assertSeqSize<args, 8> (); try { locateFunctor(my_classifier(a1, a2, a3, a4, a5, a6, a7, a8))(a1, a2, a3, a4, a5, a6, a7, a8); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4, a5, a6, a7, a8); }
    }
    
    void operator() (arg1_type a1, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
    {
      assertSeqSize<args, 9> (); try { locateFunctor(my_classifier(a1, a2, a3, a4, a5, a6, a7, a8, a9))(a1, a2, a3, a4, a5, a6, a7, a8, a9); } catch(const std::exception & e) { my_default_fn(a1, a2, a3, a4, a5, a6, a7, a8, a9); }
    }
    
  //@}
  
private:

  FunctorType       my_default_fn;
  FunctorMap        my_fns;
  ClassifierFunctor my_classifier;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
