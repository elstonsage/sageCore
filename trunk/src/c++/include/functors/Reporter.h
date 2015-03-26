#ifndef FUNCTORS_REPORTER_H
#define FUNCTORS_REPORTER_H

#include <string>
#include <iostream>
#include <sstream>
#include "boost/type_traits/add_reference.hpp"
#include "boost/mpl/push_front.hpp"
#include "functors/Misc.h"
#include "functors/ArgList.h"
#include "functors/is_void.h"
#include "functors/const_ref.h"
#include "functors/make_function.h"

namespace SAGE     {
namespace FUNCTORS {

namespace Reporter_Private
{

extern int indent;

//=========================================================================================
//
//                   DefaultPreMessageReporter
//
//=========================================================================================

///
/// The DefaultPreMessageReporter is used by the Reporter by default. It only reports the name of the executed functor.
template<typename ARGS>
class DefaultPreMessageReporter
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type  NameType;

  std::string operator() (NameType name)                                                                                                                 const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2)                                                                                                   const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3)                                                                                     const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4)                                                                       const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5)                                                         const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6)                                           const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7)                             const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8)               const { return "Beginning " + name + "-> "; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const { return "Beginning " + name + "-> "; }
};

//=========================================================================================
//
//                   DefaultPostMessageReporter - return type IS void
//
//=========================================================================================

template<bool RETURN_IS_VOID, typename ARGS> class DefaultPostMessageReporter;

template<typename ARGS>
class DefaultPostMessageReporter<true, ARGS> // Return IS void
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type NameType;
  
  std::string operator() (NameType name)                                                                                                                 const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2)                                                                                                   const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3)                                                                                     const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4)                                                                       const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5)                                                         const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6)                                           const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7)                             const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8)               const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const { return "finished " + name; }
};

//=========================================================================================
//
//                   DefaultPostMessageReporter - return type NOT void
//
//=========================================================================================

template<typename ARGS>
class DefaultPostMessageReporter<false, ARGS> // Return is NON void
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type ResultType;
  typedef arg2_type NameType;

  std::string operator() (ResultType r, NameType name)                                                                                                                 const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg3_type a3)                                                                                                   const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4)                                                                                     const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5)                                                                       const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6)                                                         const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7)                                           const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8)               const { return "finished " + name; }
  std::string operator() (ResultType r, NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const { return "finished " + name; }
};


//=========================================================================================
//
//                   CompletePreMessageReporter
//
//=========================================================================================

template<typename ARGS>
class CompletePreMessageReporter
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type          NameType;

  std::string operator() (NameType name) const 
  { 
    std::ostringstream s; s << "Beginning " + name + "--> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2) const 
  { 
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3) const 
  { 
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4) const 
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5) const
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ", " << a7 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ", " << a7 << ", " << a8 << ") --> "; return s.str(); 
  }

  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
  {
    std::ostringstream s; s << "Beginning " + name + " (" << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ", " << a7 << ", " << a8 << ", " << a9 << ") --> "; return s.str(); 
  }

};


//=========================================================================================
//
//                   CompletePostMessageReporter - return type IS void
//
//=========================================================================================

template<bool RETURN_IS_VOID, typename ARGS> class CompletePostMessageReporter;

template<typename ARGS>
class CompletePostMessageReporter<true, ARGS> // Return IS void
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type NameType;
  
  std::string operator() (NameType name)                                                                                                                 const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2)                                                                                                   const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3)                                                                                     const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4)                                                                       const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5)                                                         const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6)                                           const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7)                             const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8)               const { return "finished " + name; }
  std::string operator() (NameType name, arg2_type a2, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const { return "finished " + name; }
};

//=========================================================================================
//
//                   CompletePostMessageReporter - return type NOT void
//
//=========================================================================================

template<typename ARGS>
class CompletePostMessageReporter<false, ARGS> // Return is NOT void
{
public:

  DECLARE_ARGLIST;
  
  typedef arg1_type ResultType;
  typedef arg2_type NameType;

  std::string operator() (ResultType r, NameType name) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }
  
  std::string operator() (ResultType r, NameType name, arg3_type a3) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }

  std::string operator() (ResultType r, NameType name, arg3_type a3, arg4_type a4, arg5_type a5, arg6_type a6, arg7_type a7, arg8_type a8, arg9_type a9) const
  {
    std::ostringstream s; s << "finished with value " << r; return s.str();
  }
};

//=========================================================================================
//
//                   Executor - return type IS void
//
//=========================================================================================

template<bool RETURN_IS_VOID> class Executor { };

template<> class Executor<true> // Return IS void
{
public:

  template<typename REPORTER> static void go(const REPORTER & reporter)
  {
    assertSeqSize<typename REPORTER::args, 0> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName());

    reporter.getFunctor() ();
      
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName()) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1)
  {
    assertSeqSize<typename REPORTER::args, 1> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1);

    reporter.getFunctor() (a1);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2)
  {
    assertSeqSize<typename REPORTER::args, 2> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2);

    reporter.getFunctor() (a1, a2);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3)
  {
    assertSeqSize<typename REPORTER::args, 3> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3);

    reporter.getFunctor() (a1, a2, a3);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3) << std::endl;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4)
  {
    assertSeqSize<typename REPORTER::args, 4> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4);

    reporter.getFunctor() (a1, a2, a3, a4);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4) << std::endl;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5)
  {
    assertSeqSize<typename REPORTER::args, 5> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5);

    reporter.getFunctor() (a1, a2, a3, a4, a5);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6)
  {
    assertSeqSize<typename REPORTER::args, 6> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6);

    reporter.getFunctor() (a1, a2, a3, a4, a5, a6);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7)
  {
    assertSeqSize<typename REPORTER::args, 7> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7);

    reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7, typename REPORTER::arg8_type a8)
  {
    assertSeqSize<typename REPORTER::args, 8> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8);

    reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7, a8);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8) << std::endl;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7, typename REPORTER::arg8_type a8,
                                                                        typename REPORTER::arg9_type a9)
  {
    assertSeqSize<typename REPORTER::args, 9> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8, a9);

    reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7, a8, a9);
    
    reporter.getOs() << reporter.getPostMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8, a9) << std::endl;
  }
};

//=========================================================================================
//
//                   Executor - return type NOT void
//
//=========================================================================================


template<> class Executor<false> // Return NOT void
{
public:

  template<typename REPORTER> static void go(const REPORTER & reporter)
  {
    assertSeqSize<typename REPORTER::args, 0> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName());
  
    typename REPORTER::result_type r(reporter.getFunctor()());
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName()) << std::endl;
    
    return r;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1)
  {
    assertSeqSize<typename REPORTER::args, 1> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1);

    typename REPORTER::result_type r(reporter.getFunctor() (a1));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1) << std::endl;
    
    return r;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2)
  {
    assertSeqSize<typename REPORTER::args, 2> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2) << std::endl;
    
    return r;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3)
  {
    assertSeqSize<typename REPORTER::args, 3> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3) << std::endl;
    
    return r;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4)
  {
    assertSeqSize<typename REPORTER::args, 4> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4) << std::endl;
    
    return r;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5)
  {
    assertSeqSize<typename REPORTER::args, 5> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4, a5));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4, a5) << std::endl;
    
    return r;
  }

  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6)
  {
    assertSeqSize<typename REPORTER::args, 6> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4, a5, a6));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4, a5, a6) << std::endl;

    return r;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7)
  {
    assertSeqSize<typename REPORTER::args, 7> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4, a5, a6, a7) << std::endl;
    
    return r;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7, typename REPORTER::arg8_type a8)
  {
    assertSeqSize<typename REPORTER::args, 8> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7, a8));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8) << std::endl;

    return r;
  }
  
  template<typename REPORTER> static void go(const REPORTER & reporter, typename REPORTER::arg1_type a1, typename REPORTER::arg2_type a2,
                                                                        typename REPORTER::arg3_type a3, typename REPORTER::arg4_type a4,
                                                                        typename REPORTER::arg5_type a5, typename REPORTER::arg6_type a6,
                                                                        typename REPORTER::arg7_type a7, typename REPORTER::arg8_type a8,
                                                                        typename REPORTER::arg9_type a9)
  {
    assertSeqSize<typename REPORTER::args, 9> (); 

    reporter.getOs() << reporter.getPreMessageFunctor()(reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8, a9);

    typename REPORTER::result_type r(reporter.getFunctor() (a1, a2, a3, a4, a5, a6, a7, a8, a9));
    
    reporter.getOs() << reporter.getPostMessageFunctor()(r, reporter.getName(), a1, a2, a3, a4, a5, a6, a7, a8, a9) << std::endl;
    
    return r;
  }
};


} // End namespace Reporter_Private

//=========================================================================================
//
//                   Reporter
//  
//=========================================================================================

template<typename R, typename ARGS = boost::mpl::vector<> >
class Reporter
{
private:

  typedef Reporter_Private::Executor<is_void<R>::value> ExecutorType;

public:

  /// @name Miscellaneous
  //@{
  
    DECLARE_ARGLIST;

    typedef typename make_function<R, ARGS>::type FunctorType;
    typedef typename FunctorType::result_type     result_type;
    
  //@}

  /// @name Pre- and Post- message argument lists
  //@{
  
    typedef typename boost::mpl::push_front<ARGS, const std::string &>::type   PreMessageArgs; // name + args

    // If the return type is void then the post message arglist is the same as the pre message list
    // If the return type is non-void then we have to prepend a const reference to the return type to the post message arglist

    typedef typename boost::mpl::if_<typename is_void<R>::type,
              PreMessageArgs,
              typename boost::mpl::push_front<PreMessageArgs, typename const_ref<result_type>::type>::type>::type PostMessageArgs;
              
  //@}
  
  /// @name Pre- and Post- message functors
  //@{

    typedef typename make_function<std::string, PreMessageArgs>::type  PreMessageFunctor;
    typedef typename make_function<std::string, PostMessageArgs>::type PostMessageFunctor; 
    
  //@}

  /// @name Default pre- and post- message functors
  //@{

    typedef Reporter_Private::DefaultPreMessageReporter  <PreMessageArgs>                   DefaultPreMessageReporter;
    typedef Reporter_Private::DefaultPostMessageReporter <is_void<R>::value, PreMessageArgs> DefaultPostMessageReporter;
    
  //@}

  /// @name Constructors / operators
  //@{

    Reporter() :
      my_name ("[no name]"),
      my_os   (&std::cout),
      my_pre  (DefaultPreMessageReporter  ()),
      my_post (DefaultPostMessageReporter ())
    { }

    ///
    /// Default constuctor.
    /// This constructor will cause the Reporter to use a DefaultEvaluator (every step is assumed to be successful).
    Reporter(const FunctorType & fn, const std::string & name = "[no name]") : 
      my_name (name), 
      my_os   (&std::cout),
      my_fn   (fn), 
      my_pre  (DefaultPreMessageReporter  ()),
      my_post (DefaultPostMessageReporter ())
    { }

    ///
    /// Copy constructor.
    Reporter(const Reporter & other) : 
      my_name (other.my_name),
      my_os   (other.my_os),
      my_fn   (other.my_fn),
      my_pre  (other.my_pre),
      my_post (other.my_post) 
    { }

    ///
    /// Assignment operator.
    Reporter& operator=(const Reporter & other) 
    { 
      my_name = other.my_name; 
      my_os   = other.my_os;
      my_fn   = other.my_fn; 
      my_pre  = other.my_pre;
      my_post = other.my_post;
      
      return *this; 
    }
    
  //@}
  
  /// @name Setting components
  //@{
  
    ///
    /// Sets the Reporter to report all information. This requires that all function arguments be
    /// streamable to an ostream.
    void setReportAll() 
    { 
      my_pre  = Reporter_Private::CompletePreMessageReporter  <PreMessageArgs>                    ();
      my_post = Reporter_Private::CompletePostMessageReporter <is_void<R>::value, PostMessageArgs> ();
    }
    
    ///
    /// Sets the Reporter to report only the name of the executed functor.
    void setReportBasic() 
    { 
      my_pre  = DefaultPreMessageReporter  ();
      my_post = DefaultPostMessageReporter ();
    }

    ///
    /// Sets the name of the reported functor.
    void setName(const std::string & name) { my_name = name; }
    
    ///
    /// Sets the stream to which reporting information should be directed.
    void setStream(std::ostream & os) { my_os = os; }

    ///
    /// Sets the functor to report / execute.
    void setFunctor(const FunctorType & fn) { my_fn = fn; }

    ///
    /// Sets the premessage functor.
    void setPreMessageReporter(const PreMessageFunctor & fn) { my_pre = fn; }
  
    ///
    /// Sets the postmessage functor.
    void setPostMessageReporter(const PostMessageFunctor & fn) { my_post = fn; }
  
  
  //@}
  
  /// @name Accessing components
  //@{
  
    ///
    /// Returns the name of the reported functor.
    const std::string & getName() const { return my_name; }
    
    ///
    /// Returns the ostream to which messages will be reported.
    std::ostream & getOs() const { return *my_os; }
    
    ///
    /// Returns a const reference to the reported functor.
    const FunctorType & getFunctor() const { return my_fn; }
    
    ///
    /// Returns a const reference to the premessage functor.
    const PreMessageFunctor & getPreMessageFunctor() const { return my_pre; }
    
    ///
    /// Returns a const reference to the postmessage functor.
    const PostMessageFunctor & getPostMessageFunctor() const { return my_post; }
    
  //@}

  /// @name Invoking the sequence
  //@{

    result_type operator() ()                                         const { return ExecutorType::go(*this);             }  
    result_type operator() (arg1_type a1)                             const { return ExecutorType::go(*this, a1);         }  
    result_type operator() (arg1_type a1, arg2_type a2)               const { return ExecutorType::go(*this, a1, a2);     }  
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3) const { return ExecutorType::go(*this, a1, a2, a3); }      
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4)                             const { return ExecutorType::go(*this, a1, a2, a3, a4); }      
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4, arg5_type a5)               const { return ExecutorType::go(*this, a1, a2, a3, a4, a5); }
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4, arg5_type a5, arg6_type a6) const { return ExecutorType::go(*this, a1, a2, a3, a4, a5, a6); }
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4, arg5_type a5, arg6_type a6,
                            arg7_type a7)                             const { return ExecutorType::go(*this, a1, a2, a3, a4, a5, a6, a7); }
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4, arg5_type a5, arg6_type a6,
                            arg7_type a7, arg8_type a8)               const { return ExecutorType::go(*this, a1, a2, a3, a4, a5, a6, a7, a8); }
    result_type operator() (arg1_type a1, arg2_type a2, arg3_type a3,
                            arg4_type a4, arg5_type a5, arg6_type a6,
                            arg7_type a7, arg8_type a8, arg9_type a9) const { return ExecutorType::go(*this, a1, a2, a3, a4, a5, a6, a7, a8, a9); }

  //@}
  
private:

          std::string        my_name;
  mutable std::ostream*      my_os;
          FunctorType        my_fn;
          PreMessageFunctor  my_pre;
          PostMessageFunctor my_post;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
