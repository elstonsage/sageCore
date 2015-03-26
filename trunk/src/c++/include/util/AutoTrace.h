#ifndef UTIL_AUTO_TRACE_H
#define UTIL_AUTO_TRACE_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>

#include "util/StringUtils.h"

#define NOARGS \
std::string("")

#define ARGLIST1(arg1) \
std::string(#arg1) << "=" << (arg1)

#define ARGLIST2(arg1, arg2) \
std::string(#arg1) << "=" << (arg1) << " " << \
std::string(#arg2) << "=" << (arg2)

#define ARGLIST3(arg1, arg2, arg3) \
std::string(#arg1) << "=" << (arg1) << " " << \
std::string(#arg2) << "=" << (arg2) << " " << \
std::string(#arg3) << "=" << (arg3)

#define ARGLIST4(arg1, arg2, arg3, arg4) \
std::string(#arg1) << "=" << (arg1) << " " << \
std::string(#arg2) << "=" << (arg2) << " " << \
std::string(#arg3) << "=" << (arg3) << " " << \
std::string(#arg4) << "=" << (arg4)

#define ARGLIST5(arg1, arg2, arg3, arg4, arg5) \
std::string(#arg1) << "=" << (arg1) << " " << \
std::string(#arg2) << "=" << (arg2) << " " << \
std::string(#arg3) << "=" << (arg3) << " " << \
std::string(#arg4) << "=" << (arg4) << " " << \
std::string(#arg5) << "=" << (arg5)

#ifdef AUTOTRACE_ON

#define AUTOTRACE(arglist, f) \
f \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_STATIC(classname, arglist, f) \
f \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(classname), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT1(arglist, f) \
f  \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT2(arglist, f, init2) \
f , init2 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT3(arglist, f, init2, init3) \
f , init2, init3 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT4(arglist, f, init2, init3, init4) \
f , init2, init3, init4 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT5(arglist, f, init2, init3, init4, init5) \
f , init2, init3, init4, init5 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT6(arglist, f, init2, init3, init4, init5, init6) \
f , init2, init3, init4, init5, init6 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#define AUTOTRACE_CT7(arglist, f, init2, init3, init4, init5, init6, init7) \
f , init2, init3, init4, init5, init6, init7 \
std::ostringstream __at_arglist; \
__at_arglist << arglist; \
SAGE::UTIL::AutoTrace __auto_trace__(typeid(*this), #f, __FILE__, __LINE__, __at_arglist.str());

#else

#define AUTOTRACE(arglist, f) f

#define AUTOTRACE_STATIC(classname, arglist, f) f

#define AUTOTRACE_CT1(arglist, f)                                           f
#define AUTOTRACE_CT2(arglist, f, init2)                                    f, init2
#define AUTOTRACE_CT3(arglist, f, init2, init3)                             f, init2, init3
#define AUTOTRACE_CT4(arglist, f, init2, init3, init4)                      f, init2, init3, init4
#define AUTOTRACE_CT5(arglist, f, init2, init3, init4, init5)               f, init2, init3, init4, init5
#define AUTOTRACE_CT6(arglist, f, init2, init3, init4, init5, init6)        f, init2, init3, init4, init5, init6
#define AUTOTRACE_CT7(arglist, f, init2, init3, init4, init5, init6, init7) f, init2, init3, init4, init5, init6, init7

#endif

namespace SAGE {
namespace UTIL {

/// \brief Assists in debugging code at runtime
///
/// \par Introduction
///
/// Welcome to the Wild Word of AutoTrace! Although AutoTrace's internal implementation is very complex,
/// making use of AutoTrace is relatively simple. But first, let's review what the purpose of AutoTrace is...
///
/// First, though, let me clarify what I mean by 'AutoTrace'. As a developer, you will never use the AutoTrace
/// class directly. Instead, there is a series of preprocessor macros available for using the AutoTrace features.
///
/// \par Purpose
///
/// The AutoTrace system makes it possible to track member function entry & exit, function arguments, and data
/// members during runtime. For any given member function that has been properly encapsulated in an AutoTrace
/// macro, at runtime all executions of the function can be reported to stdout.
///
/// \par Getting started
///
/// In order to use AutoTrace, you can do it the easy way, or the slightly harder way:
///
/// The easy way:
///
/// (1) Include "util/AutoTrace.h" in your program.
///
/// (2) Add the AUTOTRACE macro to the definition of any class function that you want to track.
///
/// (3) Run make with 'AUTOTRACE=ON'. SAGE's make system has been configured to process this option.
///
/// The (slightly harder) way:
///
/// (1) Include "util/AutoTrace.h" in your program.
///
/// (2) Add -DAUTOTRACE_ON to your compiler's flags.
///
/// (3) Add the AUTOTRACE macro to the definition of any class function that you want to track.
///
/// (4) Link the libutil.a library to your program.
///
/// \par Adding the AUTOTRACE macro(s) in the simplest form
///
/// Consider the following code snippet:
///
/// \code
/// class MyClass
/// {
/// public:
///   int addOne(int y);
/// };
///
/// int MyClass::addOne(int y)
/// {
///   return y + 1;
/// }
///
/// int main()
/// {
///   MyClass m;
///   std::cout << m.addOne(3);
/// }
/// \endcode
///
/// When compiled and executed, you'll get the following output:
/// \verbatim
/// 4
/// \endverbatim
///
/// You can make this function AutoTrace-able by adding the AUTOTRACE macro as such:
///
/// \code
/// class MyClass
/// {
/// public:
///   int addOne(int y);
/// };
///
/// int AUTOTRACE(NOARGS, MyClass::addOne(int y)
/// {)
///   return y + 1;
/// }
///
/// int main()
/// {
///   MyClass m;
///   std::cout << m.addOne(3);
/// }
/// \endcode
///
/// When compiled with the -DAUTOTRACE_ON flag and executed, you'll get the following output:
/// \verbatim
///   [Begin] MyClass::addOne (foo.cpp, line XX)
///   [End]   MyClass::addOne (foo.cpp, line XX)
/// 4
/// \endverbatim
///      
/// Please note that the AUTOTRACE macro begins \b after the return type of the function
/// defition, and \b ends immeidately after the open brace of the function body.
///
/// AUTOTRACE takes two arguments: The first is information about which function arguments
/// should be reported, and the second is the portion of the function definition from the class
/// name to the first opening brace.
///
/// In the above example, for simplicity's sake, the first argument to AUTOTRACE has been
/// passed as NOARGS. This simply tells AutoTrace not to report any information about
/// the arguments passed to the function.
///
/// \par Printing function arguments with AutoTrace
///
/// AutoTrace can also print values of function arguments at runtime. This is accomplished by
/// the first argument to the AUTOTRACE macro. In the above example, the first argument was 
/// NOARGS; this tells AutoTrace not to report any information about function arguments
/// at runtime.
///
/// If you want to print function arguments, however, you can do so with the ARGLISTx macros, where
/// 'x' is the number of arguments to print. Let's consider our original example:
///
/// \code
/// int AUTOTRACE(NOARGS, MyClass::addOne(int y)
/// {)
///   return y + 1;
/// }
/// \endcode
///
/// In this example, at runtime you'll see when MyClass::addOne() is executed. But what if you want
/// to see that value of 'y', the argument passed to addOne() ?  This is accomlished in the following
/// snippet:
///
/// \code
/// int AUTOTRACE(ARGLIST1(y), MyClass::addOne(int y)
/// {)
///   return y + 1;
/// }
/// \endcode
///
/// When this is executed, the runtime information will now include information regarding the value
/// of 'y':
///
/// \verbatim
///   [Begin] MyClass::addOne y=3 (foo.cpp, line XX)
///   [End]   MyClass::addOne y=3 (foo.cpp, line XX)
/// 4
/// \endverbatim
///
/// Note that the first argument to the AUTOTRACE macro has been changed to ARGLIST1(y). Since we only
/// want to print one argument, we use ARGLIST1 (note the '1' at the end of the macro). Then we pass
/// one argument to ARGLIST1 (the name of the argument to print).
///
/// Please note that only variables that can be directed to an outputstream can be printed. If for some parameter
/// 'foo' there is no operator<< defined that can properly direct a foo instance to an outputstream you
/// will get a compile-time error.
///
/// \par What if I also want to print out data members along with function arguments?
///
/// No problem! AutoTrace doesn't distinguish between function arguments from data members. As long as the
/// named variable in the ARGLISTx macro is available to the function body normally, it can be printed
/// with AutoTrace.
///
/// Let's amend the original example to illustrate this feature. We'll add a data member 'my_x' to MyClass, and
/// initialize it in MyClass's constructor to 1.
///
/// \code
/// class MyClass
/// {
/// public:
///   MyClass() { my_x = 1; }
///   int addOne(int y);
/// private:
///   int my_x;
/// };
///
/// int AUTOTRACE(ARGLIST1(y), MyClass::addOne(int y)
/// {)
///   return y + 1;
/// }
///
/// int main()
/// {
///   MyClass m;
///   std::cout << m.addOne(3);
/// }
/// \endcode
///
/// When compiled with the -DAUTOTRACE_ON flag and executed, you'll get the following output:
/// \verbatim
///   [Begin] MyClass::addOne y=3 (foo.cpp, line XX)
///   [End]   MyClass::addOne y=3 (foo.cpp, line XX)
/// 4
/// \endverbatim
///
/// Now let's modify the AUTOTRACE macro so that we can also get the value of 'my_x' whenever addOne()
/// is called:
///
/// \code
/// int AUTOTRACE(ARGLIST2(my_x, y), MyClass::addOne(int y)
/// {)
///   return y + 1;
/// }
/// \endcode
///
/// Now we'll get the values for 'my_x' and 'y' at runtime:
/// \verbatim
///   [Begin] MyClass::addOne my_x=1 y=3 (foo.cpp, line XX)
///   [End]   MyClass::addOne my_x=1 y=3 (foo.cpp, line XX)
/// 4
/// \endverbatim
///
/// \par Indentation in runtime output
///
/// You may be wondering if there's some way to progressively indent the runtime output as functions
/// call other functions (that call other functions, etc.). AutoTrace automatically increases the runtime
/// indent by two characters every time an AUTOTRACE-ed function enters, and decrements the indent
/// when that function exits. Take a look at src/c++/util/test_autotrace.cpp to get an example that demonstrates
/// this feature. You'll see that at runtime, test_autotrace produces the following output:
///
/// \verbatim
///  A::go, y=4 (test_autotrace.cpp, line 35) [Begin]
///    C::C (test_autotrace.cpp, line 23) [Begin]
///    C::C (test_autotrace.cpp, line 23) [End]
///  A::go, y=4 (test_autotrace.cpp, line 35) [End]
/// \endverbatim
///
/// \par Constructors
///
/// There is a special macro for constructors with initializers in them. It's called AUTOTRACE_CTx, where x is the
/// number of initializers in the constructor.
///
/// Consider the following code:
/// \code
/// Foo::Foo() : my_x(0), my_y(1), my_z(3) { /* ... */ }
/// \endcode
///
/// The normal AUTOTRACE macros won't work for it! Instead, you should use AUTOTRACE_CT3 (because there are 3
/// initializers). Here's what your code should look like:
///
/// \code
/// AUTOTRACE_CT3(Foo::Foo() : my_x(0), my_y(1), my_z(3) {) /* ... */ }
/// \endcode
///
/// Now it will work!
///
/// \par Static functions
///
/// Static member functions will NOT compile properly with the normal AUTOTRACE macro. Instead, use the
/// AUTOTRACE_STATIC macro:
///
/// AUTOTRACE_STATIC(classname, arglist, function definition)
///
/// For example:
/// \code
/// class C
/// {
/// public:
///   static void foo();
/// };
///    
/// void
/// AUTOTRACE_STATIC(C, NOARGS,
/// C::foo() 
/// {)
///   /* ... */
/// }
/// \endcode
///
/// \par Class Masks
///
/// Ok, so let's say you've now incorporated AutoTrace into a whole lot of code. Let's say you've incorporated
/// it, in fact, into your entire program. Now whenever you run your program, you get a LOT of AutoTrace output.
/// So much output that it has become unreadable. You want to filter out certain classes from the output.
/// Is there a way?
///
/// Yes! You can filter out classes at compile with the static functions AutoTrace::excludeClass() and
/// AutoTrace::includeClass(). Each function takes class type. 
/// At runtime, AutoTrace will include / exclude classes based on your use of those functions.
///
/// If you use includeClass(), that implies that all classes should be EXCLUDED BY DEFAULT, and that only
/// the named classes in includeClass() should be included in the output.
///
/// If you use excludeClass(), that implies that all classes should be INCLUDED BY DEAFULT, and that only
/// the named classes in excludeClass() should be excluded from output.
///
/// You can ONLY use includeClass() or excludeClass() in your program, since the implicit assumptions of
/// their behavior contradict each other. If you try to use both of them in your program, AutoTrace will
/// abort execution.
///
/// Let's look at an example:
/// \code
/// class C { public: C(); }
/// class B { public: B(); }
/// class A { public: A(); }
///
/// AUTOTRACE(NOARGS, C::C() {) }
/// AUTOTRACE(NOARGS, B::B() {) C c; }
/// AUTOTRACE(NOARGS, A::A() {) B b; }
///
/// int main() { A a; }
/// \endcode
///
/// The above code, when run, will produce output something like this:
/// \verbatim
/// A::A() ... BEGIN
///   B::B() ... BEGIN
///     C::C() ... BEGIN
///     C::C() ... END
///   B::B() ... END
/// A::A() ... END
/// \endverbatim
///
/// Now let's say we only want to see output for no classes EXCEPT class B. We'll change the main() function:
/// \code
/// int main() { SAGE::UTIL::AutoTrace::includeClass<B>(); A a; }
/// \endcode
///
/// Now the program's output will look like this:
/// \verbatim
/// B::B() ... BEGIN
/// B::B() ... END
/// \endverbatim
///
/// Ok, now let's say we only want to see output for all classes EXCEPT class B. We'll change the main() function:
/// \code
/// int main() { SAGE::UTIL::AutoTrace::excludeClass<B>(); A a; }
/// \endcode
///
/// Now the program's output will look like this:
/// \verbatim
/// A::A() ... BEGIN
///   C::C() ... BEGIN
///   C::C() ... END
/// A::A() ... END
/// \endverbatim
///
class AutoTrace
{
public:

  /// \internal
  /// Traces certain key information about the calling function.
  ///
  /// \par The 'code' parameter
  ///
  /// Given some function:
  /// \code
  /// int someclass::foo(int x) { return x + 2; }
  /// \endcode
  ///
  /// The 'code' portion to be extracted and sent to this constructor is:
  /// \verbatim
  /// someclass::foo(int x) {
  /// \endverbatim
  ///
  /// Note that the 'code' portion BEGINS with the class name, and ENDS with the open bracket for the function
  /// body.
  ///
  /// \param file The name of the file of the calling function
  /// \param line The line number of the calling function in its respective file
  /// \param arglist A string containing information about arguments passed to the calling function
  AutoTrace(const std::type_info & tid, const std::string & code, const std::string & file, int line, const std::string & arglist);
  
  /// \internal
  /// Destructor.
  ~AutoTrace();
  
  /// @name Runtime class masks
  //@{
  
    ///
    /// See the section on class masks above.
    ///
    /// Template parameter CLASS_TYPE : The type of class to include
    ///
    template<class CLASS_TYPE> static bool includeClass()
    {
      if(classes_to_exclude.size())
      {
        std::cout << "AUTOTRACE: Error -- Cannot explicitly include a class, since all classes have been implicitly included." << std::endl;
        exit(0);
      }

      classes_to_include.push_back(typeid(CLASS_TYPE).name());

      std::cout << "AUTOTRACE: Including " << typeid(CLASS_TYPE).name() << std::endl;

      return true;
    }

    ///
    /// See the section on class masks above.
    ///
    /// Template parameter CLASS_TYPE : The type of class to exclude
    ///
    template<class CLASS_TYPE> static bool excludeClass()
    {
      if(classes_to_include.size())
      {
        std::cout << "AUTOTRACE: Error -- Cannot explicitly exclude a class, since all classes have been implicitly excluded." << std::endl;
        exit(0);
      }

      std::cout << "AUTOTRACE: Excluding " << typeid(CLASS_TYPE).name() << std::endl;

      classes_to_exclude.push_back(typeid(CLASS_TYPE).name());

      return true;
    }
    
  //@}

private:

  ///
  /// Renders the contents of this AutoTrace as a string (with no line breaks!)
  std::string toString() const;

  ///
  /// Verifies that this class is in the classlist mask
  bool isInClassList();

  std::string my_classname;
  std::string my_funcname;
  std::string my_file;
  std::string my_arglist;
  int my_line;

  static int indent;

  static std::vector<std::string> classes_to_include;
  static std::vector<std::string> classes_to_exclude;


};

inline AutoTrace::AutoTrace(const std::type_info & tid, const std::string & code, const std::string & file, int line, const std::string & arglist) :
  my_file    (file),
  my_arglist (arglist),
  my_line    (line)
{
  my_classname = tid.name();

  std::string composite_name = code.substr(0, code.find('('));

  size_t last_colon_idx = composite_name.find_last_of(":");
  
  my_funcname = composite_name.substr(last_colon_idx + 1);

  if(!isInClassList())
    return;

  std::cout << std::setw(indent * 2) << " " << toString() << " [Begin]" << std::endl;

  indent++;
}

inline AutoTrace::~AutoTrace()
{
  if(!isInClassList())
    return;

  indent--;
    
  std::cout << std::setw(indent * 2) << " " << toString() << " [End]" << std::endl;
}

inline std::string 
AutoTrace::toString() const
{
  std::ostringstream s;
    
  s << my_classname << "::" << my_funcname << (my_arglist != "" ? ", " + my_arglist  : "") << " (" << my_file << ", line " << my_line << ")";
  
  return s.str();
}


inline bool 
AutoTrace::isInClassList()
{
  if(classes_to_include.size())
  {
    for(size_t i = 0; i < classes_to_include.size(); ++i)
      if(my_classname == classes_to_include[i])
        return true;
        
    return false;
  }
  else if(classes_to_exclude.size())
  {
    for(size_t i = 0; i < classes_to_exclude.size(); ++i)
      if(my_classname == classes_to_exclude[i])
        return false;
        
    return true;
  }
  else
  {
    return true;
  }
}

} // End namespace UTIL
} // End namespace SAGE

#endif
