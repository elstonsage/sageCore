#ifndef CONTAINERS_ANYVECTOR_H
#define CONTAINERS_ANYVECTOR_H

#include "boost/type_traits/is_same.hpp"
#include "boost/any.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/mpl/contains.hpp"
#include "boost/operators.hpp"
#include "boost/iterator/filter_iterator.hpp"
#include "boost/static_assert.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>

namespace SAGE {

/// \brief Helper structs for AnyVector
namespace AnyVector_Private {

typedef std::vector<boost::any> vector_type;

struct NoType {};

template<typename T> struct isT { bool operator() (const boost::any & a) const { return a.type() == typeid(T); } };

} // End namespace AnyVector_Private

class AnyVectorItrs
{
public:

  /// @name Type-specific Iterators
  //@{
  
    template<typename T> class ConstIterator;
    
    template<typename T> class Iterator : public boost::bidirectional_iterator_helper<Iterator<T>, T, T *, T &>
    {
      friend class ConstIterator<T>;
      
      public:
      
        Iterator(AnyVector_Private::vector_type::iterator begin, AnyVector_Private::vector_type::iterator end) { my_itr = filter_itr(begin, end); }

        Iterator              ()                             {                                              }
        Iterator              (const Iterator & other)       { my_itr = other.my_itr;                       }
        Iterator & operator=  (const Iterator & other)       { my_itr = other.my_itr; return *this;         }
        T        & operator*  ()                       const { return *(boost::any_cast<T>(&*my_itr));      }
        T        * operator-> ()                       const { return   boost::any_cast<T>(&*my_itr);       }
        Iterator & operator++ ()                             { ++my_itr; return *this;                      }
        Iterator & operator-- ()                             { --my_itr; return *this;                      }
        bool       operator== (const Iterator & other) const { return my_itr == other.my_itr;               }
        bool       operator<  (const Iterator & other) const { return my_itr <  other.my_itr;               }
 
      private:
      
        typedef boost::filter_iterator<AnyVector_Private::isT<T>, AnyVector_Private::vector_type::iterator> filter_itr;

        filter_itr my_itr;
    };
  
    template<typename T> class ConstIterator : public boost::bidirectional_iterator_helper<ConstIterator<T>, T, const T *, const T &>
    {
      public:
      
        ConstIterator(AnyVector_Private::vector_type::const_iterator begin, AnyVector_Private::vector_type::const_iterator end) { my_itr = filter_itr(begin, end); }

        ConstIterator              ()                                  {                                              }
        ConstIterator              (const ConstIterator & other)       { my_itr = other.my_itr;                       }
        ConstIterator              (const Iterator<T>   & other)       { my_itr = other.my_itr;                       }
        ConstIterator & operator=  (const ConstIterator & other)       { my_itr = other.my_itr; return *this;         }
        const T       & operator*  ()                            const { return *(boost::any_cast<T>(&*my_itr));      }
        const T       * operator-> ()                            const { return   boost::any_cast<T>(&*my_itr);       }
        ConstIterator & operator++ ()                                  { ++my_itr; return *this;                      }
        ConstIterator & operator-- ()                                  { --my_itr; return *this;                      }
        bool            operator== (const ConstIterator & other) const { return my_itr == other.my_itr;               }
        bool            operator<  (const ConstIterator & other) const { return my_itr <  other.my_itr;               }
 
      private:
      
        typedef boost::filter_iterator<AnyVector_Private::isT<T>, AnyVector_Private::vector_type::const_iterator> filter_itr;

        filter_itr my_itr;
    };
  
  //@}
};

/// \brief A vector-like container that allows up a variable number of types
///
/// \par Introduction
///
/// An AnyVector works just like a normal vector, except you can insert any kind
/// of object you want. By default, it allows you to insert any object. You can
/// optionally specify a list of allowable types. With this option, the AnyVector
/// will fail at compile-time if user attempts to insert a disallowed object.
/// Furthermore, you can iterate across the contents of the AnyVector based on which type you want.
///
/// \par A quick example
///
/// \code
///
/// AnyVector<> a;
/// a.push_back(2);
/// a.push_back(2.5);
/// a.push_back(std::string("hello!"));
///
/// AnyVector<>::
///
/// AnyVector<int> b;
/// b.push_back(2.5); // Compiler fails - double not allowed in AnyVector<int>.
///
/// \endcode
///
/// \par How to restrict insertable types
///
/// You specify which types are allowed to be inserted by listing them as template arguments
/// at the point of instantiation. If you don't specify any template parameters (empty angle brackets,
/// that is), the AnyVector will allow anything to be inserted. By listing the allowable types, however,
/// the AnyVector will restrict which types can be inserted.
///
/// For example:
/// \code AnyVector<int, float, std::string> a; \endcode
/// The above code creates an AnyVector that allows only int's, float's, and string's to be inserted.
/// Keep in mind that any attempt to insert a disallowed type will fail at COMPILE TIME, not RUNTIME.
///
/// \par Iterating across contents
///
/// You can iterate across the contents of the vector on a type-by-type basis. The iterator, that is, lets
/// you iterate across a subset of the vector contents. The Iterator (or ConstIterator) is templatized 
/// on the type to which it will dereference. For example:
///
/// \code
///
/// AnyVector<int, float> a; a.push_back(1); a.push_back(3); a.push_back(1.5); a.push_back(3.5);
///
/// for(AnyVectorItrs::Iterator<int> itr = a.begin<int>(); itr != a.end<int>(); ++a) std::cout << *itr;
/// for(AnyVectorItrs::Iterator<double> itr = a.begin<double>(); itr != a.end<double>(); ++a) std::cout << *itr;
///
/// \endcode
///
/// The above code will print out "131.53.5" at runtime.
///
template<typename T0 = AnyVector_Private::NoType,
         typename T1 = AnyVector_Private::NoType,
         typename T2 = AnyVector_Private::NoType,
         typename T3 = AnyVector_Private::NoType,
         typename T4 = AnyVector_Private::NoType,
         typename T5 = AnyVector_Private::NoType,
         typename T6 = AnyVector_Private::NoType,
         typename T7 = AnyVector_Private::NoType,
         typename T8 = AnyVector_Private::NoType,
         typename T9 = AnyVector_Private::NoType>
         
class AnyVector
{
private:

    typedef AnyVector_Private::vector_type vector_type;

    typedef boost::mpl::vector<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> typelist;

    static const bool typelist_is_empty =
      boost::is_same<
        typelist, 
        boost::mpl::vector<AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType,
                           AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType, AnyVector_Private::NoType>
        >::value;
    
public:

  /// @name Inserting objects in the vector
  //@{

    ///
    /// Add the T instance to the back of the vector.
    template<typename T> void push_back(const T & t) {  assertObjAllowed<T> (); my_anys.push_back(t); }

  //@}
    
  /// @name Type-based iteration
  //@{
  
    ///
    /// Returns a non-const begin iterator across type T in the vector.
    template<typename T> AnyVectorItrs::Iterator<T> begin() { return AnyVectorItrs::Iterator<T>(my_anys.begin(), my_anys.end()); }    

    ///
    /// Returns a non-const end iterator across type T in the vector.
    template<typename T> AnyVectorItrs::Iterator<T> end() { return AnyVectorItrs::Iterator<T>(my_anys.end(), my_anys.end()); }
    
    ///
    /// Returns a const begin iterator across type T in the vector.
    template<typename T> AnyVectorItrs::ConstIterator<T> begin() const { return AnyVectorItrs::ConstIterator<T>(my_anys.begin(), my_anys.end()); }
    
    ///
    /// Returns a const end iterator across type T in the vector.
    template<typename T> AnyVectorItrs::ConstIterator<T> end() const {  return AnyVectorItrs::ConstIterator<T>(my_anys.end(), my_anys.end()); }
    
  //@}
  
  /// @name Type-based counts
  //@{
  
    ///
    /// Retruns the number of T instances in the vector.
    template<typename T> size_t count() const
    {
      size_t cnt = 0;

      for(size_t i = 0; i < my_anys.size(); ++i)
        if(boost::any_cast<T>(&(my_anys[i])))
          ++cnt;
          
      return cnt;
    }
  
  //@}

  /// @name Single-element extraction
  //@{
  
    ///
    /// Returns whether or not this objct has one and only one instance of OBJ_TYPE.
    ///
    template<typename T> bool hasOne() const
    {
      return count<T>() == 1;
    }
  
    ///
    /// Returns a non-const reference to the one OBJ_TYPE instance.
    ///
    /// NOTE: Throws an exception if object is not present. Use hasOneObj() to check first.
    ///
    template<typename T> T & getOnly()
    {
      if(hasOne<T>())
        return *begin<T>();
      else
        throw std::exception();
    }
   
    ///
    /// Returns a const reference to the one OBJ_TYPE instance.
    ///
    /// NOTE: Throws an exception if object is not present. Use hasOneObj() to check first.
    ///
    template<typename T> const T & getOnly() const
    {
      if(hasOne<T>())
        return *begin<T>();
      else
        throw std::exception();
    }
   
  //@}

  /// @name Misc
  //@{
  
    ///
    /// Indicates whether or not a child element can be cast as the given type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// Template parameter OBJ_TYPE: The type to check for
    /// \param index The absolute index number of the child element (in absolute insert order)
    template<typename T> bool isType(size_t index) const
    {
      if(index >= size())
      {
        throw std::exception();
      }
      else 
      {
        return my_anys[index].type() == typeid(T);
      }
    }
    
    ///
    /// Returns the absolute i'th object cast as type T.
    template<typename T> T & getAbs(size_t i) { return *boost::any_cast<T>(&(my_anys[i])); }

    ///
    /// Returns the absolute i'th object cast as type T.
    template<typename T> const T & getAbs(size_t i) const { return *boost::any_cast<T>(&(my_anys[i])); }
    
    ///
    /// Returns the i'th T added to the vector.
    template<typename T> T & getRel(size_t i) { AnyVectorItrs::Iterator<T> itr = begin<T>(); for(size_t j = 0; j < i; ++j, ++itr) { } return *itr; }

    ///
    /// Returns the i'th T added to the vector.
    template<typename T> const T & getRel(size_t i) const { AnyVectorItrs::ConstIterator<T> itr = begin<T>(); for(size_t j = 0; j < i; ++j, ++itr) { } return *itr; }

  //@}
  
  /// @name Proxy interface to the vector
  //@{
  
    ///
    /// Returns the i'th boost::any from the vector (non-const).
    boost::any & operator[] (size_t i) { return my_anys[i]; }

    ///
    /// Returns the i'th boost::any from the vector (const).
    const boost::any & operator[] (size_t i) const { return my_anys[i]; }

    ///
    /// Returns a non-const begin iterator across the wrappered std::vector<boost::any>.
    vector_type::iterator begin() { return my_anys.begin(); }

    ///
    /// Returns a non-const end iterator across the wrappered std::vector<boost::any>.
    vector_type::iterator end() { return my_anys.end(); }
  
    ///
    /// Returns a const begin iterator across the wrappered std::vector<boost::any>.
    vector_type::const_iterator begin() const { return my_anys.begin(); }

    ///
    /// Returns a const end iterator across the wrappered std::vector<boost::any>.
    vector_type::const_iterator end() const { return my_anys.end(); }
    
    ///
    /// Returns the total size of the wrappered std::vector<boost::any>.
    size_t size() const { return my_anys.size(); }
    
    ///
    /// Clears the vector.
    void clear() { my_anys.clear(); }
  
  //@}
  
  /// @name Direct access to the vector
  //@{
  
    vector_type & getVector() { return my_anys; }
    
    const vector_type & getVector() const { return my_anys; }
  
  
  //@}

private:

  ///
  /// Statically asserts (compile-time) that the T is allowed to be inserted).
  template<typename T> void assertObjAllowed() const
  {
    typedef typename boost::mpl::contains<typelist, T>::type T_is_in_list;
      
    typedef typename boost::mpl::or_<T_is_in_list, typename boost::mpl::bool_<typelist_is_empty>::type >::type proceed;
      
    BOOST_STATIC_ASSERT(( proceed::value )) ;
  }
    
  vector_type my_anys;
};

} // End namespace SAGE

#endif
