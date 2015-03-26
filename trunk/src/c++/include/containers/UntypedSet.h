#ifndef UNTYPED_SET_H
#define UNTYPED_SET_H

#include "boost/bind.hpp"
#include "boost/static_assert.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/operators.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/mpl/contains.hpp"
#include "boost/mpl/count.hpp"
#include "boost/mpl/eval_if.hpp"
#include "boost/mpl/for_each.hpp"
#include "boost/mpl/find_if.hpp"
#include "boost/type_traits/is_same.hpp"
#include "boost/type_traits/is_convertible.hpp"
#include "boost/mpl/or.hpp"
#include "util/TypeInfo.h"
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

namespace SAGE {

/// \brief Stores a set of objects of any type(s)
///
/// \par Introduction
///
/// The basic idea of the UntypedSet is that it lets you store
/// an arbitrary set of data. By arbitrary, I mean that the objects stored
/// in the UntypedSet can be of *any* type. Having placed the objects
/// in the map, you can then access then by type-specific object 
/// iterators.
///
/// \par Example
///
/// For instance, let's say you've got a bunch of Car objects and a bunch of Bike
/// objects. You'd like to store them all in a single container. The common solution
/// is to derive both Car and Bike from a base class such as Vehicle with a bunch
/// of virtual functions, and then store 
/// a vector of Vehicle pointers. Then you can loop across your vector, invoke the virtual
/// functions on each pointer, and the correctly derived version will be invoked.
///
/// It's a working solution, but not too satisfactory. What if you have a lot of unrelated
/// objects that you need to store together? What if, instead of Car and Bike (slightly related
/// objects), you have Car's, Bottle's, MessageFactory's, President's, and std::vector<std::pair<Foo, Bar> >'s ?
/// It would be seriously difficult to force them all to be derived from a common base object. And even
/// if they were, you couldn't impose a standard virtual function table on all of them.
///
/// Instead, you can use the UntypedSet! You can place as many objects of different types
/// in the UntypedSet, and extract them type-by-type later on. Consider the following code:
///
/// \code
///
/// UntypedSet m;
/// class foo { public: void doFoo(); }; 
/// class bar { public: void doBar(); };
/// m.insert(foo()); 
/// m.insert(bar()); 
/// m.insert(foo());
/// for(UntypedSet::Iterator<foo> i = m.begin<foo>(); i != m.end<foo>(); ++i)
///   i->doFoo();
/// for(UntypedSet::Iterator<bar> i = m.begin<bar>(); i != m.end<bar>(); ++i)
///   i->doBar();
///
/// \endcode
///
/// In the above example, we insert two instances of 'foo' and one instance of 'bar' into the
/// UntypedSet. Then we iterate across all foo's, invoking doFoo() on each, and across
/// all bar's, invoking doBar() on each.
///
/// \par More on iteration
///
/// Iteration is done via Iterator<...>'s. The Iterator is templatized on what kind of
/// object you want to iterate across. When you dereference it, you get a reference to that
/// object type.
///
/// For instance, if you want to iterate across all objects of type 'thingy', you enter:
/// \code
/// UntypedSet m;
/// for(UntypedSet::Iterator<thingy> i = m.begin<thingy>(); i != m.end<thingy>(); ++i)
///   ;
/// \endcode
///
/// \par Sorting
///
/// Objects of a particular type can be sorted according to user-given criterion. By default, they
/// are sorted by the order in which they are inserted in to the container. Please see addSorter()
/// for more information.
class UntypedSet
{
public:

  typedef boost::shared_ptr<void> ObjShPtr;

  typedef std::pair<int, ObjShPtr> PairedObj;

  //===============
  // Functor stuff:
  //===============

  class FunctorInternalBase
  { public: virtual bool operator() (const void *, const void *) const = 0; };

  typedef boost::shared_ptr<FunctorInternalBase> FIBShPtr;

  template<class FUNCTOR_TYPE, class INPUT_TYPE>
  class FunctorInternal : public FunctorInternalBase
  {
    public:

      FunctorInternal            (const FUNCTOR_TYPE   & functor) : my_functor(functor) {}
      FunctorInternal            (const FunctorInternal & other) { my_functor = other.my_functor; }
      FunctorInternal& operator= (const FunctorInternal & other) { my_functor = other.my_functor; return *this; }

      virtual bool operator() (const void * t1, const void * t2) const
      {
        return my_functor(*(const INPUT_TYPE*)t1, *(const INPUT_TYPE*)t2);
      }

    private:
      FUNCTOR_TYPE my_functor;
  };

  class FunctorWrapper
  {
  public:

    class DefaultFunctor : public FunctorInternalBase
    { public: virtual bool operator() (const void * t1, const void * t2) const { return t1 < t2; } };

    explicit FunctorWrapper(bool sort_by_insert_order = false)
    { 
      my_sort = sort_by_insert_order;
      my_f    = FIBShPtr(new DefaultFunctor()); 
    }
    
    FunctorWrapper            (const FunctorWrapper & other) { my_f = other.my_f; my_sort = other.my_sort; }
    FunctorWrapper& operator= (const FunctorWrapper & other) { my_f = other.my_f; my_sort = other.my_sort; return *this; }

    explicit FunctorWrapper(FIBShPtr f) : my_f(f), my_sort(false)  {}

    bool operator() (const PairedObj & t1, const PairedObj & t2) const
    {
      return my_sort ? t1.first < t2.first : my_f->operator()(t1.second.get(), t2.second.get());
    }

  private:
    FIBShPtr my_f;
    bool     my_sort;
  };


  typedef std::set<PairedObj, FunctorWrapper> PairedObjSet;

  //===========================================
  // Iterator & ConstIterator - across ObjSet's
  //===========================================

  template<typename OBJ_TYPE> class ConstIterator;

  template<typename OBJ_TYPE>
  struct Iterator : public boost::bidirectional_iterator_helper<Iterator<OBJ_TYPE>, OBJ_TYPE, OBJ_TYPE *, OBJ_TYPE &>
  {
    friend class UntypedSet;
    friend class ConstIterator<OBJ_TYPE>;

  public:
    Iterator              ()                             {                                              }
    Iterator              (const Iterator & other)       { my_itr = other.my_itr;                       }
    Iterator & operator=  (const Iterator & other)       { my_itr = other.my_itr; return *this;         }
    OBJ_TYPE & operator*  ()                       const { return *((OBJ_TYPE*)(my_itr->second.get())); }
    OBJ_TYPE * operator-> ()                       const { return  ((OBJ_TYPE*)(my_itr->second.get())); }
    Iterator & operator++ ()                             { ++my_itr; return *this;                      }
    Iterator & operator-- ()                             { --my_itr; return *this;                      }
    bool       operator== (const Iterator & other) const { return my_itr == other.my_itr;               }
    bool       operator<  (const Iterator & other) const { return my_itr <  other.my_itr;               }

  private:

    explicit Iterator(PairedObjSet::iterator itr) { my_itr = itr; }
    PairedObjSet::iterator my_itr;
  };

  template<typename OBJ_TYPE>
  struct ConstIterator : public boost::bidirectional_iterator_helper<ConstIterator<OBJ_TYPE>, OBJ_TYPE, const OBJ_TYPE *, const OBJ_TYPE &>
  {
    friend class UntypedSet;

  public:
    ConstIterator               ()                                       {                                              }
    ConstIterator               (const ConstIterator      & other)       { my_itr = other.my_itr;                       }
    ConstIterator               (const Iterator<OBJ_TYPE> & other)       { my_itr = other.my_itr;                       }
    ConstIterator  & operator=  (const ConstIterator      & other)       { my_itr = other.my_itr; return *this;         }
    const OBJ_TYPE & operator*  ()                                 const { return *((OBJ_TYPE*)(my_itr->second.get())); }
    const OBJ_TYPE * operator-> ()                                 const { return  ((OBJ_TYPE*)(my_itr->second.get())); }
    ConstIterator  & operator++ ()                                       { ++my_itr; return *this;                      }
    ConstIterator  & operator-- ()                                       { --my_itr; return *this;                      }
    bool             operator== (const ConstIterator      & other) const { return my_itr == other.my_itr;               }
    bool             operator<  (const ConstIterator      & other) const { return my_itr <  other.my_itr;               }

  private:
 
    explicit ConstIterator(PairedObjSet::const_iterator itr) { my_itr = itr; }
    PairedObjSet::const_iterator my_itr;
  };

  //==============
  // DataMap type
  //==============

  typedef std::map<UTIL::TypeInfo, PairedObjSet> DataMap;

  typedef std::vector<std::pair<UTIL::TypeInfo, ObjShPtr> > ObjVector;

  //============
  // AbsIterator
  //============

  template<typename OBJ_TYPE> class AbsConstIterator;

  template<typename OBJ_TYPE>
  struct AbsIterator : public boost::bidirectional_iterator_helper<UntypedSet::AbsIterator<OBJ_TYPE>, OBJ_TYPE, OBJ_TYPE *, OBJ_TYPE &>
  {
    friend class UntypedSet;
    friend class AbsConstIterator<OBJ_TYPE>;

  public:
    AbsIterator              ()                                {                                      }
    AbsIterator              (const AbsIterator & other)       { my_itr = other.my_itr;               }
    AbsIterator & operator=  (const AbsIterator & other)       { my_itr = other.my_itr; return *this; }
    OBJ_TYPE    & operator*  ()                          const { return *(OBJ_TYPE*)(my_itr->second.get());  }
    OBJ_TYPE    * operator-> ()                          const { return  (OBJ_TYPE*)(my_itr->second.get());  }
    bool          operator== (const AbsIterator & other) const { return my_itr == other.my_itr;       }
    bool          operator<  (const AbsIterator & other) const { return my_itr < other.my_itr;        }
    AbsIterator & operator++ ()                                { ++my_itr; return *this;              }
    AbsIterator & operator-- ()                                { --my_itr; return *this;              }

  private:

    explicit AbsIterator(ObjVector::iterator itr) : my_itr(itr) {}

    ObjVector::iterator my_itr;
  };

  //=================
  // AbsConstIterator
  //=================

  template<typename OBJ_TYPE>
  struct AbsConstIterator : public boost::bidirectional_iterator_helper<UntypedSet::AbsConstIterator<OBJ_TYPE>, OBJ_TYPE, OBJ_TYPE *, OBJ_TYPE &>
  {
    friend class UntypedSet;

  public:
    AbsConstIterator              ()                                          {                                      }
    AbsConstIterator              (const AbsConstIterator      & other)       { my_itr = other.my_itr;               }
    AbsConstIterator              (const AbsIterator<OBJ_TYPE> & other)       { my_itr = other.my_itr;               }
    AbsConstIterator & operator=  (const AbsConstIterator      & other)       { my_itr = other.my_itr; return *this; }
    OBJ_TYPE         & operator*  ()                                    const { return *(OBJ_TYPE*)(my_itr->second.get());  }
    OBJ_TYPE         * operator-> ()                                    const { return  (OBJ_TYPE*)(my_itr->second.get());  }
    bool               operator== (const AbsConstIterator      & other) const { return my_itr == other.my_itr;       }
    bool               operator<  (const AbsConstIterator      & other) const { return my_itr < other.my_itr;        }
    AbsConstIterator & operator++ ()                                          { ++my_itr; return *this;              }
    AbsConstIterator & operator-- ()                                          { --my_itr; return *this;              }

  private:

    explicit AbsConstIterator(ObjVector::const_iterator itr) : my_itr(itr) {}

    ObjVector::const_iterator my_itr;
  };

  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    UntypedSet()
    {
      clear();
    }

    ///
    /// Copy constructor.
    UntypedSet(const UntypedSet & other)
    {
      my_max_counts = other.my_max_counts;
      my_data       = other.my_data;
      my_error_set  = other.my_error_set;
      my_objs  = other.my_objs;
    }

    ///
    /// Assignment operator.
    UntypedSet& operator=(const UntypedSet & other)
    {
      my_max_counts = other.my_max_counts;
      my_data       = other.my_data;
      my_error_set  = other.my_error_set;
      my_objs  = other.my_objs;
      
      return *this;
    }

  //@}

  /// @name Object insertion
  //@{

    ///
    /// Sets the maximum number of OBJ_TYPE's allowed to be inserted in this object.
    ///
    template<typename OBJ_TYPE> void setMaxCount(int count)
    {
      my_max_counts[typeid(OBJ_TYPE)] = count;
    }
    
    ///
    /// Returns the maximum number of OBJ_TYPE's allowed to be inserted in this object.
    ///
    template<typename OBJ_TYPE> int getMaxCount() const
    {
      std::map<UTIL::TypeInfo, int>::const_iterator i = my_max_counts.find(typeid(OBJ_TYPE));
      
      if(i == my_max_counts.end())
        return std::numeric_limits<int>::infinity();
    }

    ///
    /// Indicates that objects of type OBJ_TYPE are to be sorted according to a specific
    /// sorting functor.
    ///
    /// The functor is required to have the following public function defined:
    ///
    ///    bool operator() (const OBJ_TYPE & o1, const OBJ_TYPE & o2) const;
    ///
    /// It should return true if o1 is less than o2, or false otherwise.
    ///
    /// Please note: This function can only be used if no objects of type OBJ_TYPE have 
    /// been added yet. If it is invoked \b after the insertion of such objects, an
    /// exception will be thrown.
    ///
    template<typename OBJ_TYPE, typename SORTER_TYPE>
    void addSorter(const SORTER_TYPE & sorter)
    {
      DataMap::iterator set_itr = my_data.find(typeid(OBJ_TYPE));
    
      if(set_itr == my_data.end()) // If no objects have been added yet, we can add a sorter:
      {
        FunctorWrapper functor(FIBShPtr(new FunctorInternal<SORTER_TYPE, OBJ_TYPE>(sorter)));

        my_data.insert(std::make_pair(UTIL::TypeInfo::create<OBJ_TYPE>(), PairedObjSet(functor)));
      }
      else
      {
        throw std::exception();
      }
    }

    ///
    /// Inserts the object.
    ///
    /// Note: If the maximum object count is exceed, this function will simply not insert the given object.
    ///
    template<typename OBJ_TYPE> void insert(const OBJ_TYPE & obj)
    {
      // Check the maximum allowable count first:
      
      std::map<UTIL::TypeInfo, size_t>::const_iterator i = my_max_counts.find(typeid(OBJ_TYPE));

      if(i == my_max_counts.end() || ((i != my_max_counts.end()) && (count<OBJ_TYPE>() < i->second)))
      {
        // If the set hasn't been created yet, create it with the insert-order-sorter:
        
        if(my_data.find(typeid(OBJ_TYPE)) == my_data.end())
          sortByInsertOrder<OBJ_TYPE> ();
      
        // Create the shared pointer and stick it in the set:
        
        ObjShPtr sh_ptr(new OBJ_TYPE(obj));
        
        my_data[typeid(OBJ_TYPE)].insert(PairedObj(count<OBJ_TYPE>(), sh_ptr));
        
        my_objs.push_back(std::make_pair(UTIL::TypeInfo::create<OBJ_TYPE>(), sh_ptr));
      }
    }
    
  //@}

  /// @name Single-child extraction
  //@{
  
    ///
    /// Returns whether or not this objct has one and only one instance of OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> bool hasOneObj() const
    {
      return count<OBJ_TYPE>() == 1;
    }
  
    ///
    /// Returns a non-const reference to the one OBJ_TYPE instance.
    ///
    /// NOTE: Throws an exception if object is not present. Use hasOneObj() to check first.
    ///
    template<typename OBJ_TYPE> OBJ_TYPE& getOnlyObj()
    {
      if(hasOneObj<OBJ_TYPE>())
        return *begin<OBJ_TYPE>();
      else
        throw std::exception();
    }

    ///
    /// Returns a const reference to the one OBJ_TYPE instance.
    ///
    /// NOTE: Throws an exception if object is not present. Use hasOneObj() to check first.
    ///
    template<typename OBJ_TYPE> const OBJ_TYPE& getOnlyObj() const
    {
      if(hasOneObj<OBJ_TYPE>())
        return *begin<OBJ_TYPE>();
      else
        throw std::exception();
    }

  //@}

  /// @name Indexed extraction
  //@{
  
    ///
    /// Indicates whether or not a child element can be cast as the given type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// Template parameter OBJ_TYPE: The type to check for
    /// \param index The absolute index number of the child element (in absolute insert order)
    template<typename OBJ_TYPE> bool isOfType(size_t index) const
    {
      if(index >= size()) 
      {
        throw std::exception();
      }
      else 
      {
        return my_objs[index].first == typeid(OBJ_TYPE);
      }
    }

    ///
    /// Returns the i'th element as the requested type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// \param index The index number of the element, corresponding to the absolute order in which elements were
    /// inserted into the UntypedSet.
    ///
    template<typename OBJ_TYPE> const OBJ_TYPE & getAbs(size_t index) const
    {
      if(index >= size()) throw std::exception();
      else                return *(OBJ_TYPE*)(my_objs[index].second.get());
    }
  
    ///
    /// Returns the i'th element as the requested type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// \param index The index number of the element, corresponding to the absolute order in which elements were
    /// inserted into the UntypedSet.
    ///
    template<typename OBJ_TYPE> OBJ_TYPE & getAbs(size_t index)
    {
      if(index >= size()) throw std::exception();
      else                return *(OBJ_TYPE*)(my_objs[index].second.get());
    }
  
    ///
    /// Returns the i'th object of OBJ_TYPE (where i'th is determined by the sort order for OBJ_TYPE).
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    template<typename OBJ_TYPE> const OBJ_TYPE & get(size_t index) const
    {
      if(index >= count<OBJ_TYPE>())
      {
        throw std::exception();
      }
      else
      {
        ConstIterator<OBJ_TYPE> itr = begin<OBJ_TYPE> ();
        
        if(index > 0)
          for(size_t i = 0; i < index; ++i, ++itr) { }
          
        return *itr;
      }
    }

    ///
    /// Returns the i'th object of OBJ_TYPE (where i'th is determined by the sort order for OBJ_TYPE).
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    template<typename OBJ_TYPE> OBJ_TYPE & get(size_t index)
    {
      if(index >= count<OBJ_TYPE>())
      {
        throw std::exception();
      }
      else
      {
        Iterator<OBJ_TYPE> itr = begin<OBJ_TYPE> ();
        
        if(index > 0)
          for(size_t i = 0; i < index; ++i, ++itr) { }
          
        return *itr;
      }
    }

  //@}

  /// @name Subset iteration
  //@{
  
    ///
    /// Returns a const begin iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> ConstIterator<OBJ_TYPE> begin() const
    { 
      DataMap::const_iterator i = my_data.find(typeid(OBJ_TYPE));
      
      return ConstIterator<OBJ_TYPE>(i == my_data.end() ? my_error_set.begin() : i->second.begin());
    }
  
    ///
    /// Returns a const end iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> ConstIterator<OBJ_TYPE> end() const
    { 
      DataMap::const_iterator i = my_data.find(typeid(OBJ_TYPE));
      
      return ConstIterator<OBJ_TYPE>(i == my_data.end() ? my_error_set.end() : i->second.end());
    }
    
    ///
    /// Returns a begin iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> Iterator<OBJ_TYPE> begin()
    { 
      DataMap::iterator i = my_data.find(typeid(OBJ_TYPE));
      
      return Iterator<OBJ_TYPE>(i == my_data.end() ? my_error_set.begin() : i->second.begin());
    }
    
    ///
    /// Returns a end iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> Iterator<OBJ_TYPE> end()
    { 
      DataMap::iterator i = my_data.find(typeid(OBJ_TYPE));
      
      return Iterator<OBJ_TYPE>(i == my_data.end() ? my_error_set.end() : i->second.end());
    }

  //@}
  
  /// @name Complete set iteration
  //@{
  
    ///
    /// Returns a non-const begin iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> AbsIterator<OBJ_TYPE> beginAbs() { return AbsIterator<OBJ_TYPE> (my_objs.begin()); }
  
    ///
    /// Returns a non-const end iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> AbsIterator<OBJ_TYPE> endAbs() { return AbsIterator<OBJ_TYPE> (my_objs.end()); }
    
    ///
    /// Returns a const begin iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> AbsConstIterator<OBJ_TYPE> beginAbs() const { return AbsConstIterator<OBJ_TYPE> (my_objs.begin()); }
  
    ///
    /// Returns a const end iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> AbsConstIterator<OBJ_TYPE> endAbs() const
    { 
      return AbsConstIterator<OBJ_TYPE> (my_objs.end());
    }
    
  //@}
  
  /// @name Miscellaneous
  //@{

    ///
    /// Clears the entire data structure.
    void clear()
    {
      my_data       . clear();
      my_max_counts . clear();
      my_objs  . clear();
    }

    ///
    /// Returns the total number of elements in the map.
    DataMap::size_type size() const
    {
      return my_objs.size();
    }

    ///
    /// Returns the number of elements of type OBJ_TYPE in the map.
    template<typename OBJ_TYPE> DataMap::size_type count() const
    { 
      DataMap::const_iterator i = my_data.find(typeid(OBJ_TYPE));
      
      return i == my_data.end() ? 0 : i->second.size();
    }
    
  //@}

private:

    ///
    /// Indicates that objects of type OBJ_TYPE should be sorted according to the order in which
    /// they are inserted.
    ///
    template<typename OBJ_TYPE>
    void sortByInsertOrder()
    {
      DataMap::iterator set_itr = my_data.find(typeid(OBJ_TYPE));
    
      if(set_itr == my_data.end()) // If no objects have been added yet, we can add a sorter:
      {
        my_data.insert(std::make_pair(UTIL::TypeInfo::create<OBJ_TYPE>(), PairedObjSet(FunctorWrapper(true))));
      }
      else
      {
        throw std::exception();
      }
    }



  DataMap                          my_data;
  PairedObjSet                     my_error_set;
  std::map<UTIL::TypeInfo, size_t> my_max_counts;
  ObjVector                        my_objs;
};

/// \internal
/// \brief Helper struct for doing implicit type conversion on insertable types
///
/// Here's the deal: If you've got a TypedSet that only allows std::string's, and
/// someone tries to insert "foo", the compiler may interpret "foo" as a char array,
/// NOT as a std::string. This is a problem!
///
/// The ConvertAndInsert struct, however, handles it. Here's what happens: The TypedSet's insert()
/// function checks through the allowable typelist to see if it can convert the given
/// type to any of the types in the list. If it can't, it'll invoke a compile-time
/// error (no surprise). If it can, and if the convertible type is DIFFERENT from
/// the object's type, it will create a ConvertAndInsert templatized on the
/// type it found (the 'TO_TYPE'), and invoke's the ObjInserter's convertAndInsert()
/// function. convertAndInsert() will construct a TO_TYPE out of the FROM_TYPE and insert
/// *that* into the underlying UntypedSet.
///
/// If the convertible type is the SAME as the object type, the TypedSet will create
/// a Insert struct to simply insert the object.
template<typename TO_TYPE>
struct ConvertAndInsert
{
  template<typename FROM_TYPE, typename CONTAINER_TYPE> 
  static void go(const FROM_TYPE & original_obj, CONTAINER_TYPE & container)
  {
    container.insert(TO_TYPE(original_obj));
  }
};

/// \internal
/// \brief Helper struct for inserting objects in a TypedSet
///
/// See documentation for ConvertAndInsert for more information.
template<typename TO_TYPE>
struct Insert
{
  template<typename CONTAINER_TYPE> 
  static void go(const TO_TYPE & original_obj, CONTAINER_TYPE & container)
  {
    container.insert(original_obj);
  }
};



/// \internal Helper for the UntypedSet
class DummyClass { private: DummyClass(); };

template<class T0,
         class T1 = DummyClass,
         class T2 = DummyClass,
         class T3 = DummyClass,
         class T4 = DummyClass,
         class T5 = DummyClass,
         class T6 = DummyClass,
         class T7 = DummyClass,
         class T8 = DummyClass,
         class T9 = DummyClass>

class TypedSet
{
public:

  /// \internal
  /// Boost vector of allowable insert types
  typedef boost::mpl::vector<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> allowable_types;

  /// @name Constructors / operators
  //@{
  
    ///
    /// Constructor.
    TypedSet()
    {
      my_set.clear();
    }
    
    ///
    /// Copy constructor.
    TypedSet(const TypedSet & other)
    {
      my_set = other.my_set;
    }
    
    ///
    /// Assignment operator.
    TypedSet& operator= (const TypedSet & other)
    {
      my_set = other.my_set;
      
      return *this;
    }
  
  //@}

  /// @name Inserting data
  //@{
  
    ///
    /// Sets the maximum number of OBJ_TYPE's allowed to be inserted in this object.
    template<typename OBJ_TYPE> void setMaxCount(int count)
    {
      my_set.setMaxCount<OBJ_TYPE>(count);
    }
     
    template<typename OBJ_TYPE> int getMaxCount() const
    {
      return my_set.getMaxCount<OBJ_TYPE>();
    }

    ///
    /// SORTER_TYPE must have:
    ///    bool operator< (const OBJ_TYPE &, const OBJ_TYPE &) const;
    ///
    template<typename OBJ_TYPE, typename SORTER_TYPE>
    void addSorter(const SORTER_TYPE & sorter)
    {
      my_set.addSorter<OBJ_TYPE>(sorter);
    }

    ///
    /// Inserts the object.
    ///
    template<typename OBJ_TYPE> void insert(const OBJ_TYPE & obj)
    { 
      using namespace boost;

      // Locate the position of a convertible type, if any:
      typedef typename mpl::find_if<allowable_types, is_convertible<OBJ_TYPE, mpl::_1> >::type convertible_type_pos;

      // Make sure the convertible type itr is NOT the end itr:
      BOOST_STATIC_ASSERT(( mpl::not_< is_same< convertible_type_pos, typename mpl::end<allowable_types>::type > >::value ));

      // Dereference the convertible type iterator:
      typedef typename mpl::deref<convertible_type_pos>::type convertible_type;
      
      // If the convertible type is the same as OBJ_TYPE, then just insert it. Otherwise, convert-and-insert it:
      mpl::if_<is_same<convertible_type, OBJ_TYPE>, Insert<convertible_type>, ConvertAndInsert<convertible_type> >::type::go(obj, my_set);
    }
    
    ///
    /// Inserts a series of objects into the multimap.
    /// \param begin A begin iterator
    /// \param end An end iterator
    template<typename OBJ_TYPE, typename ITERATOR_TYPE> void insert(ITERATOR_TYPE begin, ITERATOR_TYPE end)
    {
      std::for_each(begin, end, boost::bind(&insert, boost::ref(*this), _1));
    }
    
  //@}
  
  /// @name Single-child extraction
  //@{
  
    ///
    /// Please set UntypedSet reference.
    template<typename OBJ_TYPE> bool hasOneObj() const
    {
      return my_set.hasOneObj<OBJ_TYPE>();
    }
  
    ///
    /// Please set UntypedSet reference.
    template<typename OBJ_TYPE> OBJ_TYPE& getOnlyObj()
    {
      return my_set.getOnlyObj<OBJ_TYPE>();
    }

    ///
    /// Please set UntypedSet reference.
    template<typename OBJ_TYPE> const OBJ_TYPE& getOnlyObj() const
    {
      return my_set.getOnlyObj<OBJ_TYPE>();
    }

  //@}

  /// @name Indexed extraction
  //@{
  
    ///
    /// Indicates whether or not a child element can be cast as the given type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// Template parameter OBJ_TYPE: The type to check for
    /// \param index The absolute index number of the child element (in absolute insert order)
    template<typename OBJ_TYPE> bool isOfType(size_t index) const
    {
      return my_set.isOfType<OBJ_TYPE>(index);
    }
    
    ///
    /// Returns the i'th element as the requested type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// \param index The index number of the element, corresponding to the absolute order in which elements were
    /// inserted into the UntypedSet.
    ///
    template<typename OBJ_TYPE> const OBJ_TYPE & getAbs(size_t index) const
    {
      return my_set.getAbs<OBJ_TYPE>(index);
    }
  
    ///
    /// Returns the i'th element as the requested type.
    ///
    /// NOTE: Throws an exception if the index is out of range.
    ///
    /// \param index The index number of the element, corresponding to the absolute order in which elements were
    /// inserted into the UntypedSet.
    ///
    template<typename OBJ_TYPE> OBJ_TYPE & getAbs(size_t index)
    {
      return my_set.getAbs<OBJ_TYPE>(index);
    }
  
    ///
    /// See UntypedSet::get().
    ///
    template<typename OBJ_TYPE> const OBJ_TYPE & get(size_t index) const
    {
      return my_set.get<OBJ_TYPE>(index);
    }

    ///
    /// See UntypedSet::get().
    ///
    template<typename OBJ_TYPE> OBJ_TYPE & get(size_t index)
    {
      return my_set.get<OBJ_TYPE>(index);
    }

  //@}

  /// @name Subset iteration
  //@{

    ///
    /// Returns a const begin iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> UntypedSet::ConstIterator<OBJ_TYPE> begin() const
    { 
      isObjAllowed<OBJ_TYPE> ();
      
      return my_set.begin<OBJ_TYPE> ();
    }

    ///
    /// Returns a const end iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> UntypedSet::ConstIterator<OBJ_TYPE> end() const
    { 
      isObjAllowed<OBJ_TYPE> ();
      
      return my_set.end<OBJ_TYPE> ();
    }
    
    ///
    /// Returns a begin iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> UntypedSet::Iterator<OBJ_TYPE> begin()
    { 
      isObjAllowed<OBJ_TYPE> ();
      
      return my_set.begin<OBJ_TYPE> ();
    }
    
    ///
    /// Returns a end iterator for all objects of type OBJ_TYPE.
    ///
    template<typename OBJ_TYPE> UntypedSet::Iterator<OBJ_TYPE> end()
    { 
      isObjAllowed<OBJ_TYPE> ();
      
      return my_set.end<OBJ_TYPE> ();
    }
    
  //@}

  /// @name Complete set iteration
  //@{
  
    ///
    /// Returns a non-const begin iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> UntypedSet::AbsIterator<OBJ_TYPE> beginAbs()
    {
      return my_set.beginAbs<OBJ_TYPE>();
    }
  
    ///
    /// Returns a non-const end iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> UntypedSet::AbsIterator<OBJ_TYPE> endAbs()
    { 
      return my_set.endAbs<OBJ_TYPE>();
    }
    
    ///
    /// Returns a const begin iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> UntypedSet::AbsConstIterator<OBJ_TYPE> beginAbs() const
    {
      return my_set.beginAbs<OBJ_TYPE>();
    }
  
    ///
    /// Returns a const end iterator for all objects in the set.
    ///
    template<typename OBJ_TYPE> UntypedSet::AbsConstIterator<OBJ_TYPE> endAbs() const
    { 
      return my_set.endAbs<OBJ_TYPE>();
    }
    
  //@}
  
  /// @name Miscellaneous
  //@{

    ///
    /// Clears the entire data structure.
    void clear()
    {
      my_set.clear();
    }

    ///
    /// Returns the total number of elements in the map.
    UntypedSet::DataMap::size_type size() const
    {
      return my_set.size();
    }

    ///
    /// Returns the number of elements of type OBJ_TYPE in the map.
    template<typename OBJ_TYPE> UntypedSet::DataMap::size_type count() const
    { 
      return my_set.count<OBJ_TYPE> ();
    }
    
  //@}

  ///
  /// Fails at COMPILE-TIME if the child object is not allowed!
  template<typename OBJ_TYPE> static void isObjAllowed()
  {
    // Make sure OBJ_TYPE is allowed!
    /* You have tried to insert or extract a disallowed object. */ BOOST_STATIC_ASSERT((boost::mpl::contains<allowable_types, OBJ_TYPE>::value));
  }

private:

  UntypedSet my_set;
};

} /// End namespace

#endif
