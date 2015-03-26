#ifndef OUTPUT_ELEMENT_H
#define OUTPUT_ELEMENT_H

#include "error/internal_error.h"
#include "boost/shared_ptr.hpp"
#include "libxml++/libxml++.h"
#include "util/AutoTrace.h"
#include "containers/AnyVector.h"
#include "boost/mpl/if.hpp"
#include "boost/mpl/switch.hpp"
#include "boost/type_traits/remove_bounds.hpp"
#include "boost/type_traits/add_const.hpp"
#include "boost/type_traits/remove_const.hpp"
#include "boost/type_traits/is_integral.hpp"
#include "boost/type_traits/is_convertible.hpp"
#include "util/OutlineCntr.h"
#include "output/AttrMgr.h"
#include "output/Attrs.h"
#include "output/Command.h"
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace SAGE   {
namespace OUTPUT {

class RenderingRules;
class SpannedCell;
class TableColumn;
class TableRow;
class FloatingPoint;
class Section;
class Table;
class UnavailableCell;

template <typename VALUE_TYPE>
class BasicValue;

typedef BasicValue<double>      Double;
typedef BasicValue<int>         Int;
typedef BasicValue<std::string> String;

/// \internal
/// \brief Stores the information about the final derived element type, and provides insert validation functions
///
///
template<typename SELF_TYPE = AnyVector_Private::NoType>
class ElementBase
{
public:
  typedef SELF_TYPE self;
  
  /// @name Insert validation
  //@{
  
    template<typename CHILD_TYPE> bool validate(self & parent, CHILD_TYPE & child) const
    {
      return true;
    }
    
  //@}
};

///
/// This is a compiler macro for defining your own insert-validating function.
///
/// For instance, let's say you've got a 'Foo' object that can take a 'Bar' object
/// to be inserted, but the Bar object must have the same 'size' as the 'Foo' object.
/// Then your code might look
/// like this:
///
/// \code
/// ELEMENT_VALIDATOR(Foo, Bar)
/// {
///   return parent.size == child.size;
/// }
/// \endcode
///
/// Note that the parent and child elements are available as 'parent' and 'child'
/// variables.
///
#define ELEMENT_VALIDATOR(parent_type, child_type) \
template<> template<> inline bool \
ElementBase<parent_type>::validate(parent_type & parent, child_type & child) const

///
/// This is a compiler macro for defining your own command objects for elements.
///
/// Let's 
#define ELEMENT_COMMAND(element_type, command_code) \
template<> inline void \
Command<element_type, command_code>::operator() (element_type & element) const

#define ELEMENT_COMMAND1(element_type, command_code, arg_type1) \
template<> inline void \
Command1<element_type, command_code, arg_type1>::operator() (element_type & element) const

#define ELEMENT_COMMAND2(element_type, command_code, arg_type1, arg_type2) \
template<> inline void \
Command2<element_type, command_code, arg_type1, arg_type2>::operator() (element_type & element) const

/// \brief Basic interface for any output element
///
/// \par Introduction
///
/// An Element is an object with the following characteristics:
///
/// - Stores an arbitrary number of specified child types (enforced at compile-time)
///
/// - Has an arbitrary number of attributes, where attribute is a name/value pair
///
/// - Allows / disallows the "TITLE" attribute specifically (enforced at compile-time)
///
/// - Can (optionally) define its own logic for rendering its contents as text
///
/// Practically speaking, the Element class is used as a basic for creating a variety
/// of data storage components. This includes basic values (Int, Double, and String),
/// organizational components (Section), special types (RenderingRules), and
/// more complex data (Table).
///
/// \par Deriving your own component
///
/// \par Child elements
///
/// Child elements are copied into a parent element via the insert() function. This function performs
/// a complete (not shallow!) copy.
///
/// If child element need to be validated by by the parent in some way, you can write a validation function.
/// It should appear in one of your include files, and should follow the form:
///
/// \code
/// ELEMENT_VALIDATOR(parent_type, child_type)
/// {
///   // ...
///   // return true, false ?
/// }
/// \endcode
///
/// ELEMENT_VALIDATOR is a compiler macro; parent_type is the type of the parent object; child_type
/// is the type of the child to validate. Within the function, 'parent' is a non-const reference to the parent
/// object, and 'child' is a const reference to the child object. The function must return true or false (true
/// if the element may be inserted, false otherwise).
///
/// For an example, take a look at the end of Table.h, where the ELEMENT_VALIDATOR is defined for a number of
/// circumstances.
///
/// You can iterate across all child elements of a particular type via the begin() and end() functions
/// (templatized on the child element type). You can also iterate across all elements (cast as a shared base type)
/// with the beginAbs() and endAbs() functions. The beginAbs() and endAbs() functions will iterate across all
/// child elements in the absolute order in which they were insert()-ed into the parent element.
///
/// \par Attributes
///
/// \par Output formatting
///
/// You can render the contents of any Element derivation in one of two ways. If you want the contents directly,
/// you can use the toPrettyPrint(), toPlaintext(), or toXml(), each of which returns its objects contents
/// in the form in question.
///
/// Alternately, you can have the contents redirected to an output with the generatePlaintextFile(),
/// generatePrettyPrintFile(), and generateXmlFile() functions. These functions simply place the rendered contents
/// in a file.
///
template<typename SELF_TYPE = AnyVector_Private::NoType,
         typename ATTRMGR_TYPE = AttrMgr<>,
         typename T0 = AnyVector_Private::NoType, 
         typename T1 = AnyVector_Private::NoType,
         typename T2 = AnyVector_Private::NoType,
         typename T3 = AnyVector_Private::NoType,
         typename T4 = AnyVector_Private::NoType,
         typename T5 = AnyVector_Private::NoType,
         typename T6 = AnyVector_Private::NoType,
         typename T7 = AnyVector_Private::NoType,
         typename T8 = AnyVector_Private::NoType,
         typename T9 = AnyVector_Private::NoType>

class Element : public ElementBase<SELF_TYPE>, public ATTRMGR_TYPE
{
public:

  typedef SELF_TYPE self;

  typedef AnyVector<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> vector_type;

  typedef typename ATTRMGR_TYPE::AttributeMap AttributeMap;

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor
    /// \param name The typename of the derived element (such as "Table", or "Double")
    explicit Element(const std::string & name) : my_name(name)
    {
      my_vector.clear();
    }
  
    ///
    /// Copy constructor
    Element(const Element & other) : ElementBase<SELF_TYPE>(other), ATTRMGR_TYPE(other), my_vector(other.my_vector), my_name(other.my_name) { }
    
    ///
    /// Assignment operator
    Element& operator=(const Element & other)
    {
      if(this != &other)
      {
        ElementBase<SELF_TYPE>::operator=(other);
        ATTRMGR_TYPE::operator=(other);

        my_name   = other.my_name;
        my_vector = other.my_vector;
      }
      
      return *this;
    }
    
    virtual ~Element() { } 

  //@}
  
  /// @name Basic info
  //@{
  
    ///
    /// Returns the canonical name of the derived Element type. Non-word characters (colons,
    /// angle brackets) will be replaced with underscores.
    const std::string & getType() const 
    { 
      return my_name;
    }

  //@}
  
  /// @name Child element adding and iterating
  //@{

    /// \internal
    /// \brief A special struct for type-converting basic types
    ///
    /// Ok, here's the deal. We want users to be able to directly insert
    /// strings, double's, and int's WITHOUT having to wrapper them in
    /// String's, Double's, and Int's.
    ///
    /// The problem is that you can't simply have the user insert() a double,
    /// because the underlying TypedSet will cause a compile failure when
    /// it looks for 'double' in its typelist. THe typelist will have Double,
    /// but not double.
    ///
    /// How to get around this? More template meta-programming!
    ///
    /// The Element::insert() function does a compile-time check against
    /// the type of the inserted object, and choose which template specialization
    /// of addChild to create. addChild::go() then knows how to copy construct
    /// the child to be added.
    ///
    /// Woohoo!
    ///
    template<typename CHILD_TYPE>
    struct addChild
    {
      template<typename OBJ_TYPE>
      static void go(Element & parent, typename Element::vector_type & s, const OBJ_TYPE & child)
      {
        CHILD_TYPE o(child); if(parent.validate((typename Element::self &)parent, o)) s.push_back(o);
      }
    };

    ///
    /// This struct takes care of processing element commands.
    /// See SAGE::OUTPUT::Command for more information.
    template<typename PARENT_TYPE>
    struct executeCommand
    {
      template<typename CMD_TYPE>
      static void go(PARENT_TYPE & parent, typename Element::vector_type & s, const CMD_TYPE & cmd)
      {
        cmd(parent);
      }
    };

    ///
    /// Inserts the object.
    /// If you want to insert multiple objects, remember that this function returns a reference
    /// to the parent element. You can simply chain together insert statement, like so:
    ///
    /// \code
    /// TableRow r;
    /// r.insert(Double(1.0)).insert(String("Hi!"));
    /// \endcode
    ///
    /// Or with the stream operator:
    ///
    /// \code
    /// TableRow r;
    /// r << Double(1.0) << String("Hi!");
    /// \endcode
    ///
    /// \param obj The object to insert
    /// \returns A reference to the parent element.
    template<typename OBJ_TYPE> SELF_TYPE & insert(const OBJ_TYPE & obj)
    {
      using namespace boost::mpl;
      using namespace boost;

      typedef is_convertible <OBJ_TYPE, CommandBase<SELF_TYPE> > is_command;
      typedef is_same        <OBJ_TYPE, double>                  is_double;
      typedef is_integral    <OBJ_TYPE>                          treat_as_int;
      typedef is_convertible <OBJ_TYPE, std::string>             treat_as_string;

      if_<is_command, executeCommand<SELF_TYPE>,       typename
          if_<is_double, addChild<Double>,             typename
              if_<treat_as_int, addChild<Int>,         typename
                  if_<treat_as_string, addChild<String>, addChild<OBJ_TYPE> >::type 
              >::type 
          >::type 
      >::type::go((SELF_TYPE&)*this, my_vector, obj);

      return (SELF_TYPE&)*this;
    }

    ///
    /// Inserts the object (just like insert() ).
    template<typename OBJ_TYPE> self & operator<< (const OBJ_TYPE & child)
    {
      return insert(child);
    }

  //@}
  
  /// @name Underlying AnyVector
  //@{
  
    vector_type & getVector() { return my_vector; }
    
    const vector_type & getVector() const { return my_vector; }
    
  //@}

  /// @name Checking for elements with specific titles
  //@{
  
    ///
    /// Indicates whether this object has a child element whose TITLE attribute is title.
    ///
    /// Note: The OBJ_TYPE in question must have a title attribute for this function
    /// to compile.
    template<typename OBJ_TYPE> bool hasElement(const std::string & title) const
    {
      for(AnyVectorItrs::ConstIterator<OBJ_TYPE> i = my_vector.template begin<OBJ_TYPE>(); i != my_vector.template end<OBJ_TYPE>(); ++i)
        if(i->getTitle() == title)
          return true;

      return false;
    }

  //@}

private:

  // The vector_type manages the child elements
  vector_type my_vector;

  // Name of the derived type:
  std::string my_name;
};

} // End namespace OUTPUT
} // End namespace SAGE


#endif
