#ifndef UTIL_TYPEINFO_H
#define UTIL_TYPEINFO_H

#include <typeinfo>
#include <string>

namespace SAGE {
namespace UTIL {

class TypeInfo
{
  public:

    template<typename T>
    struct TypeWrapper
    {
      typedef T type;
    };

  /// @name Constructors / destructors
  //@{

    ///
    /// Constructor.
    TypeInfo() { my_info = 0; }

    ///
    /// Constructor #2.
    TypeInfo(const std::type_info & info) { my_info = &info; }

    ///
    /// Constructor #3.
    TypeInfo(const std::type_info * info) { my_info = info; }

    template<typename T> static TypeInfo create() { return TypeInfo(typeid(T)); }

    ///
    /// Copy constructor.
    TypeInfo(const TypeInfo & other) { my_info = other.my_info; }
    
  //@}

  /// @name Operators
  //@{

    TypeInfo& operator=(const TypeInfo & other) { my_info = other.my_info; return *this; }
    
    bool operator== (const TypeInfo & other) const { return  my_info && other.my_info ? my_info->operator==(*other.my_info) : false; }
    bool operator!= (const TypeInfo & other) const { return !operator==(other); }
    bool operator<  (const TypeInfo & other) const { return  before(other); }
    bool operator<= (const TypeInfo & other) const { return  before(other) || operator==(other); }
    bool operator>  (const TypeInfo & other) const { return !before(other); }
    bool operator>= (const TypeInfo & other) const { return !before(other) || operator==(other); }

  //@}

  /// Misc.
  //@{
  
    // Compatibility functions
    bool before(const TypeInfo & other) const { return my_info && other.my_info ? my_info->before(*other.my_info) : false; }

    ///
    /// Returns the name of the type, or an empty string
    /// if no type has been given yet.
    std::string getName() const { return my_info ? my_info->name() : ""; }
    
  //@}

  private:

    const std::type_info * my_info;
};


// Comparison operators

inline bool operator== (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) == a2; }
inline bool operator!= (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) != a2; }
inline bool operator<  (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) <  a2; }
inline bool operator<= (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) <= a2; }
inline bool operator>  (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) >  a2; }
inline bool operator>= (const std::type_info * a1, const TypeInfo & a2) { return TypeInfo(a1) >= a2; }

} // end namespace UTIL
} // end namespace SAGE

#endif
