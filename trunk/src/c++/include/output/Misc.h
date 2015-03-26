#ifndef OUTPUT_MISC_H
#define OUTPUT_MISC_H

#include <sys/stat.h>
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include "output/Element.h"
#include "output/RenderingRules.h"
#include <string>

namespace SAGE   {
namespace OUTPUT {

//===========
// BasicValue
//===========

/// \brief Represents a simple value
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
/// There are three specialized versions of this class: Double, Int, and String.
template<typename VALUE_TYPE>
class BasicValue : public Element<BasicValue<VALUE_TYPE>, AttrMgr<HasValue<VALUE_TYPE> > >
{
public:

  typedef Element<BasicValue<VALUE_TYPE>, AttrMgr<HasValue<VALUE_TYPE> > > Base;

  typedef VALUE_TYPE val_type;

  ///
  /// Constructor
  /// \param v The value to represent  
  explicit BasicValue(const VALUE_TYPE & v) : Base("BasicValue")
  {
    /* You have attempted to create a BasicValue other than double / int / string. */ BOOST_STATIC_ASSERT((sizeof(VALUE_TYPE) == 0));
  }

  ///
  /// Copy constructor
  BasicValue(const BasicValue & other) : Base(other) { }

  ///
  /// Assignment operator
  BasicValue& operator= (const BasicValue & other)
  {
    Base::operator=(other);
    return *this;
  }

  ///
  /// Returns the underlying represented value.
  const VALUE_TYPE & toVal() const { return Base::getValue(); }
};

/// \brief Represents a double value
typedef BasicValue<double> Double;

/// \brief Represents an int value
typedef BasicValue<int> Int;

/// \brief Represents a string value
typedef BasicValue<std::string> String;


/// Specialize the constructors for Double, Int, and String, so that the Element gets named correctly:

template<> inline BasicValue<double>      ::BasicValue(const double      & v) : Base ("Double")      { setValue(v); }
template<> inline BasicValue<int>         ::BasicValue(const int         & v) : Base ("Int")         { setValue(v); }
template<> inline BasicValue<std::string> ::BasicValue(const std::string & v) : Base ("String")      { setValue(v); }

//===========
// NamedValue
//===========

/// \brief A name paired with a value
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
/// To get the name portion of the named value, use its getTitle() function.
///
/// To get the value portion of the named value, simply cast the NamedValue as the value type.
template<class VALUE_TYPE>
class NamedValue : public Element<NamedValue<VALUE_TYPE>, AttrMgr<HasTitle, HasValue<VALUE_TYPE>, HasTooltip, HasHelpText > >
{
public:

  typedef Element<NamedValue<VALUE_TYPE>, AttrMgr<HasTitle, HasValue<VALUE_TYPE>, HasTooltip, HasHelpText > > Base;

  typedef VALUE_TYPE valtype;
  
  /// @name Constructor / operators
  //@{
  
    ///
    /// Constructor
    /// \param title The name of the NamedValue
    /// \param val The value of the NamedValue
    NamedValue(const std::string & title, const VALUE_TYPE & val) : Base("NamedValue") { Base::setTitle(title); Base::setValue(val); }

    ///
    /// Copy constructor
    NamedValue(const NamedValue  & other) : Base(other) { }

    ///
    /// Assignment operator
    NamedValue & operator= (const NamedValue & other) { Base::operator=(other); return *this; }


    ///
    /// Returns the underlying represented value.
    const VALUE_TYPE & toVal() const { return Base::getValue(); }
    
  //@}
};

/// \brief Represents a named double value.
typedef NamedValue<double> NamedDouble;

/// \brief Represents a named integer value.
typedef NamedValue<int> NamedInt;

/// \brief Represents a named string value.
typedef NamedValue<std::string> NamedString;

// Specialize the constructors for NamedDouble, NamedInt, and NamedString, so that the Element's name gets set correctly:

template<> inline NamedValue<double>      ::NamedValue(const std::string & title, const double      & val) : Base("NamedDouble") { setTitle(title); setValue(val); }
template<> inline NamedValue<int>         ::NamedValue(const std::string & title, const int         & val) : Base("NamedInt")    { setTitle(title); setValue(val); }
template<> inline NamedValue<std::string> ::NamedValue(const std::string & title, const std::string & val) : Base("NamedString") { setTitle(title); setValue(val); }

/// \brief A bulleted list
///
/// \par Introduction
///
/// Each List element represents a single bullet (as well as its (optional) sub-lists
/// and text content). Let's say you want a simple bulleted list. You'll create a 
/// List element, and stick List elements (with BulletType = BULLETED) in it. You can
/// use the static function makeBullet() to speed things up a bit. For instance:
///
/// \code
/// List l;
/// l << List::makeBullet("thing") << List::makeBullet("foo") << List::makeBullet("bar");
/// \endcode
///
/// \par Rendering
///
class List : public Element<List, AttrMgr<HasTitle, HasBulletType>, List, String>
{
public:

  typedef Element<List, AttrMgr<HasTitle, HasBulletType>, List, String> Base;
  
  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    List(const std::string & title = "", HasBulletType::BulletType e = HasBulletType::UNSPECIFIED) : Base("List")
    {
      setTitle(title);
      setBulletType(e);
    }
  
    ///
    /// Copy constructor.
    List(const List & other) : Base(other) {}
  
    ///
    /// Assignment operator.
    List& operator=(const List & other)
    {
      Base::operator=(other);
      return *this;
    }

    ///
    /// Creates a List element of bulleted type with only a title.
    static List makeBullet(const std::string & txt)
    {
      return List(txt, HasBulletType::BULLETED);
    }

  //@}
};

ELEMENT_VALIDATOR(List, List)
{
  // Make sure that all child List's have the same bullet type as the
  // first child List:
  if(parent.getVector().count<List>())
  {
    child.setBulletType(parent.getVector().begin<List>()->getBulletType());
  }
  
  return true;
}

//========
// Section
//========

/// \brief Organizes a section of related data
///
/// ALLOWABLE CHILD ELEMENTS: [Section, Table, String, NamedDouble, NamedInt, NamedString]
///
/// In order to organize sets of components, use this class. You are allowed to insert
/// Section, Table, String, NamedDouble, NamedInt, NamedString.
///
class Section : public Element<Section, AttrMgr<HasTitle, HasTooltip, HasHelpText>, Section, Table, NamedDouble, NamedInt, NamedString, List> 
{ 
public: 

  typedef Element<Section, AttrMgr<HasTitle, HasTooltip, HasHelpText>, Section, Table, NamedDouble, NamedInt, NamedString, List> Base;

  ///
  /// Constructor.
  /// \param title The title of this section
  explicit Section(const std::string & title = "") : Base ("Section") { setTitle(title); }
  
  ///
  /// Copy Constructor.
  Section(const Section & other) : Base (other) {}
  
  ///
  /// Assignment operator.
  Section& operator=(const Section & other) { Base::operator=(other); return *this; }
};


} // End namespace OUTPUT
} // End namespace SAGE

#endif
