#ifndef OUTPUT_COMMAND_H
#define OUTPUT_COMMAND_H

#include "boost/mpl/vector.hpp"
#include "boost/mpl/contains.hpp"
#include "boost/mpl/count.hpp"
#include "boost/mpl/eval_if.hpp"
#include "boost/mpl/or.hpp"
#include "boost/mpl/if.hpp"
#include "boost/mpl/switch.hpp"
#include "boost/type_traits/is_same.hpp"
#include "boost/type_traits/remove_bounds.hpp"
#include "boost/type_traits/add_const.hpp"
#include "boost/type_traits/is_convertible.hpp"

namespace SAGE   {
namespace OUTPUT {

/// \internal
/// \brief Base class for element commands
template<typename OBJ_TYPE> class CommandBase {};

/// \brief Non-specialized form for element commands
///
/// \par What is an element command?
///
/// The Command<...> class represents an action that can be sent to an Element-derived type.
/// For instance, if you look at SAGE::OUTPUT::Table, you'll see that in addition to its
/// member function SAGE::OUTPUT::Table::insertBlankRow(), there's also a typedef called
/// SAGE::OUTPUT::TABLE::INSERT_BLANK_ROW. Using this typedef, you can stream the "insert
/// blank row" command to your table:
///
/// \code
/// Table t;
/// t << (TableRow() << 5) << Table::INSERT_BLANK_ROW() << (TableRow() << "hello");
/// \endcode
///
/// \par Show me an example!
///
/// \par How can I write my own?
///

/** \brief For commands taking no arguments. */
template<typename OBJ_TYPE, 
         void (OBJ_TYPE::*MFP_NAME)() >
         
class Command0 : public CommandBase<OBJ_TYPE>
{
public:

  void operator() (OBJ_TYPE & parent) const { (parent.*MFP_NAME) (); };

};

template<typename ARG1,
         typename OBJ_TYPE, 
         void (OBJ_TYPE::*MFP_NAME)(ARG1) >
         
/** \brief For commands taking a single argument. */
class Command1 : public CommandBase<OBJ_TYPE>
{
public:

  explicit Command1(ARG1 arg1) : my_arg1(arg1) { }

  void operator() (OBJ_TYPE & parent) const { (parent.*MFP_NAME) (my_arg1); };

private:

  ARG1 my_arg1;  
};

/** \brief For commands taking two arguments. */
template<typename ARG1,
         typename ARG2,
         typename OBJ_TYPE, 
         void (OBJ_TYPE::*MFP_NAME)(ARG1, ARG2) >
         
class Command2 : public CommandBase<OBJ_TYPE>
{
public:

  Command2(ARG1 arg1, ARG2 arg2) : my_arg1(arg1) { }

  void operator() (OBJ_TYPE & parent) const { (parent.*MFP_NAME) (my_arg1, my_arg2); };

private:

  ARG1 my_arg1;  
  ARG2 my_arg2;
};


//================
// MACRO LANGUAGE
//================

#define BEGIN_COMMAND0(element_type, function_name) \
private: void _ ## function_name() {

#define BEGIN_COMMAND1(element_type, function_name, arg1) \
private: void _ ## function_name(arg1 a1) {

#define BEGIN_COMMAND2(element_type, function_name, arg1, arg2) \
private: void _ ## function_name(arg1 a1, arg2 a2) {

#define END_COMMAND0(element_type, function_name) \
} public: typedef Command0<element_type, &element_type::_ ## function_name> function_name; 

#define END_COMMAND1(element_type, function_name, arg1) \
} public: typedef Command1<arg1, element_type, &element_type::_ ## function_name> function_name; 

#define END_COMMAND2(element_type, function_name, arg1, arg2) \
} public: typedef Command2<arg1, arg2, element_type, &element_type::_ ## function_name> function_name; 

} // End namespace OUTPUT
} // End namespace SAGE

#endif
