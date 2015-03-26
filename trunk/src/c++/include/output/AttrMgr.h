#ifndef OUTPUT_ATTR_MGR_H
#define OUTPUT_ATTR_MGR_H

namespace SAGE   {
namespace OUTPUT {

/// \brief Base class for an attribute setter/getter class.
///
class AttrBase
{
public:
  
  ///
  /// Constructor.
  /// \param name The name of the attribute
  AttrBase(const std::string & name) : my_name(name), my_value("") { }
  
  ///
  /// Copy constructor.
  AttrBase(const AttrBase & other) : my_name(other.my_name), my_value(other.my_value) { }

  ///
  /// Assignment operator.
  AttrBase& operator= (const AttrBase & other) 
  { 
    if(this != &other)
    {
      my_name  = other.my_name;
      my_value = other.my_value;
    }
    
    return *this;
  }
  
  ///
  /// Returns the name of the attribute.
  const std::string & getAttrName() const { return my_name; }

  ///
  /// Returns the value of the attribute.
  const std::string & getAttrValue() const { return my_value; }

protected:
  
  ///
  /// Sets the value of the attribute.
  /// \param s The new attribute value
  void setAttrValue(const std::string & s) { my_value = s; }
  
private:

  std::string my_name;
  std::string my_value;
};

/// \brief The class for no attribute.
///
template<int i>
class NoAttr : public AttrBase
{
public:
  NoAttr() : AttrBase("") { }
  NoAttr(const NoAttr & other) : AttrBase(other) {}
  NoAttr& operator= (const NoAttr & other) { AttrBase::operator=(other); return *this; }
};

/// \brief An attribute manager
///
template<class ATTR0 = NoAttr<0>,
         class ATTR1 = NoAttr<1>,
         class ATTR2 = NoAttr<2>,
         class ATTR3 = NoAttr<3>,
         class ATTR4 = NoAttr<4>,
         class ATTR5 = NoAttr<5>,
         class ATTR6 = NoAttr<6>,
         class ATTR7 = NoAttr<7>,
         class ATTR8 = NoAttr<8>,
         class ATTR9 = NoAttr<9> >

class AttrMgr : public ATTR0, public ATTR1, public ATTR2, public ATTR3, public ATTR4, public ATTR5, public ATTR6, public ATTR7, public ATTR8, public ATTR9
{
public:

  typedef std::map<std::string, std::string> AttributeMap;

  AttrMgr() {}
  
  AttrMgr(const AttrMgr & other) : ATTR0(other), ATTR1(other), ATTR2(other), ATTR3(other), ATTR4(other), 
                                   ATTR5(other), ATTR6(other), ATTR7(other), ATTR8(other), ATTR9(other) {}
                                   
  AttrMgr& operator=(const AttrMgr & other)
  {
    ATTR0::operator=(other);
    ATTR1::operator=(other);
    ATTR2::operator=(other);
    ATTR3::operator=(other);
    ATTR4::operator=(other);
    ATTR5::operator=(other);
    ATTR6::operator=(other);
    ATTR7::operator=(other);
    ATTR8::operator=(other);
    ATTR9::operator=(other);
    
    return *this;
  }

  ///
  /// Populates an AttributeMap with name/value pairs from this AttrMgr.
  /// \param attr_map The AttributeMap to populate with values.
  void populateAttrMap(AttributeMap & attr_map) const
  {
    attr_map.clear();
    
    if(ATTR0::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR0::getAttrName(), ATTR0::getAttrValue())); }
    if(ATTR1::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR1::getAttrName(), ATTR1::getAttrValue())); }
    if(ATTR2::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR2::getAttrName(), ATTR2::getAttrValue())); }
    if(ATTR3::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR3::getAttrName(), ATTR3::getAttrValue())); }
    if(ATTR4::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR4::getAttrName(), ATTR4::getAttrValue())); }
    if(ATTR5::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR5::getAttrName(), ATTR5::getAttrValue())); }
    if(ATTR6::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR6::getAttrName(), ATTR6::getAttrValue())); }
    if(ATTR7::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR7::getAttrName(), ATTR7::getAttrValue())); }
    if(ATTR8::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR8::getAttrName(), ATTR8::getAttrValue())); }
    if(ATTR9::getAttrName() != "") { attr_map.insert(std::make_pair(ATTR9::getAttrName(), ATTR9::getAttrValue())); }
  }
};

#define CREATE_ENUM_ATTRIBUTE2(attr_name, default_val, val1, val2) \
\
class Has ## attr_name : public AttrBase                 \
{                                                           \
public:                                                     \
  enum attr_name { val1, val2 };                            \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  attr_name get##attr_name () const                 \
  {                                                         \
    return getAttrValue() == #val1 ? val1 : val2;           \
  }                                                         \
                                                            \
  void set##attr_name(attr_name a)                  \
  {                                                         \
    setAttrValue(a == val1 ? #val1 : #val2);                \
  }                                                         \
};

#define CREATE_ENUM_ATTRIBUTE3(attr_name, default_val, val1, val2, val3) \
\
class Has ## attr_name : public AttrBase                 \
{                                                           \
public:                                                     \
  enum attr_name { val1, val2, val3 };                      \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  attr_name get##attr_name () const                 \
  {                                                         \
    return getAttrValue() == #val1 ? val1 :                 \
           (getAttrValue() == #val2 ? val2 : val3);         \
  }                                                         \
                                                            \
  void set##attr_name(attr_name a)                  \
  {                                                         \
    setAttrValue(a == val1 ? #val1 :                        \
                (a == val2 ? #val2 : #val3));               \
  }                                                         \
};

#define CREATE_ENUM_ATTRIBUTE4(attr_name, default_val, val1, val2, val3, val4) \
\
class Has ## attr_name : public AttrBase                 \
{                                                           \
public:                                                     \
  enum attr_name { val1, val2, val3, val4 };                \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  attr_name get##attr_name () const                 \
  {                                                         \
    return getAttrValue() == #val1 ? val1 :                 \
           (getAttrValue() == #val2 ? val2 :                \
            (getAttrValue() == #val3 ? val3 : val4));       \
  }                                                         \
                                                            \
  void set##attr_name(attr_name a)                  \
  {                                                         \
    setAttrValue(a == val1 ? #val1 :                        \
                (a == val2 ? #val2 :                        \
                (a == val3 ? #val3 : #val4)));              \
  }                                                         \
};

#define CREATE_STRING_ATTRIBUTE(attr_name, default_val) \
\
class Has ## attr_name : public AttrBase                 \
{                                                           \
public:                                                     \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  const std::string & get##attr_name () const       \
  {                                                         \
    return getAttrValue();                                  \
  }                                                         \
                                                            \
  void set##attr_name(const std::string & a)        \
  {                                                         \
    setAttrValue(a);                                        \
  }                                                         \
};


#define CREATE_BOOL_ATTRIBUTE(attr_name, default_val) \
\
class Has ## attr_name : public AttrBase                    \
{                                                           \
public:                                                     \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  bool get##attr_name () const                      \
  {                                                         \
    return getAttrValue() == "TRUE";                        \
  }                                                         \
                                                            \
  void set##attr_name(bool a)                       \
  {                                                         \
    setAttrValue(a ? "TRUE" : "FALSE");                     \
  }                                                         \
};

#define CREATE_INT_ATTRIBUTE(attr_name, default_val) \
\
class Has ## attr_name : public AttrBase                    \
{                                                           \
public:                                                     \
                                                            \
  Has##attr_name () : AttrBase(#attr_name)                  \
  {                                                         \
    set##attr_name(default_val);                            \
  }                                                         \
                                                            \
  Has##attr_name (const Has##attr_name & other) :           \
    AttrBase(other) { }                                     \
                                                            \
  Has##attr_name & operator= (const Has##attr_name & other) \
  {                                                         \
    AttrBase::operator=(other); return *this;               \
  }                                                         \
                                                            \
  int get##attr_name () const                       \
  {                                                         \
    return atoi(getAttrValue().c_str());                    \
  }                                                         \
                                                            \
  void set##attr_name(int a)                        \
  {                                                         \
    std::ostringstream p_str; p_str << a;                   \
    setAttrValue(p_str.str());                              \
  }                                                         \
};

} // End namespace OUTPUT
} // End namespace SAGE

#endif
