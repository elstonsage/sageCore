#ifndef CONTAINERS_UNTYPED_OBJ_H
#define CONTAINERS_UNTYPED_OBJ_H

#include "util/TypeInfo.h"
#include <assert.h>
#include <iostream>

namespace SAGE {

class UntypedObj
{
public:

  /// @name Constructors
  //@{

    ///
    /// Default constructor.
    UntypedObj();
    
    UntypedObj(const UntypedObj & other);

    template<typename T> explicit UntypedObj(const T & t);
    
  //@}
  
  /// @name Assignment operators
  //@{

    UntypedObj & operator=(const UntypedObj & other);

    template<typename T> UntypedObj & operator=(const T & t);
    
  //@}
  
  /// @name Destructor
  //@{
    
    ~UntypedObj();

  //@}

  /// @name Direct access to the object
  //@{

    ///
    /// Returns \true if the underlying object is not null (true if it's valid to use)
    bool notNull() const { return my_obj; }

    ///
    /// Returns true if this object stores an object of type T; false otherwise.
    template<typename T> bool isType() const { return my_obj ? getInfo() == UTIL::TypeInfo::create<T>() : false; }

    ///
    /// Returns a non-const reference to the underlying type, cast as a T reference.
    template<typename T> T & get() { assert(my_obj); return *((T*)my_obj); }

    ///
    /// Returns a const reference to the underlying type, cast as a T reference.
    template<typename T> const T & get() const { assert(my_obj); return *((const T*)my_obj); }

  //@}  

  /// @name Checking for allowable type assignment
  //@{
  
    ///
    /// Returns a UTIL::TypeInfo object describing the type currently stored.
    /// If the underlying object hasn't yet been assigned, the behavior is undefined. You should
    /// check that my_obj is not null before invoking this function.
    UTIL::TypeInfo getInfo() const { return (this->*my_get_info_fnptr) (); }

  //@}

private:

  /// @name Assigning function pointers
  //@{

    ///
    /// Sets the obj and all function pointers to NULL.
    void reset();

    ///
    /// Sets all three function pointers to the correctly templatized versions.
    template<typename T> void assignFunctionPtrs();
    
  //@}

  /// @name Underlying templatized functions
  //@{

    ///
    /// Copies this object's contents to the other object.
    template<typename T> void copyMeToOther(UntypedObj & other) const;
    
    ///
    /// Deletes the obj if it is not null, then invokes reset().
    template<typename T> void destroyMe();

    ///
    /// Returns the TypeInfo instance for the stored type. If T is void, the behavior is undefined.
    template<typename T> UTIL::TypeInfo getInfo() const;

  //@}
  

  typedef void           (UntypedObj::*FNPTR_copyMeToOther) (UntypedObj & other) const;
  typedef void           (UntypedObj::*FNPTR_destroyMe)     ();
  typedef UTIL::TypeInfo (UntypedObj::*FNPTR_getInfo)       () const;
  
  void                * my_obj;
  FNPTR_copyMeToOther   my_copy_to_other_fnptr;
  FNPTR_destroyMe       my_destroy_me_fnptr;
  FNPTR_getInfo         my_get_info_fnptr;

};

//====================
//  INLINE FUNCTIONS
//====================

//============================================
//  CONSTRUCTOR
//============================================
inline UntypedObj::UntypedObj()
{ 
  reset();
}
    
//============================================
//  COPY CONSTRUCTOR
//============================================
inline UntypedObj::UntypedObj(const UntypedObj & other)
{
  reset();
  
  if(other.my_obj)
    (other.*(other.my_copy_to_other_fnptr))(*this);
}

//============================================
//  Constructor with OBJ
//============================================
template<typename T> 
inline UntypedObj::UntypedObj(const T & t) 
{
  my_obj = new T(t); 
  assignFunctionPtrs<T> ();
}

//============================================
//  Assignment operator with OBJ
//============================================
template<typename T> 
inline UntypedObj & UntypedObj::operator=(const T & t) 
{
  if(my_obj)
  {
    if(UTIL::TypeInfo::create<T>() == getInfo())
    {
      *((T*)my_obj) = t;

      return *this;
    }
    else
    {
      (this->*my_destroy_me_fnptr)();
    }
  }

  my_obj = new T(t);
  assignFunctionPtrs<T> ();
  
  return *this;
}

//============================================
//  Assignment operator
//============================================
inline UntypedObj & 
UntypedObj::operator=(const UntypedObj & other) 
{
  // Precheck:
  if(&other == this)
    return *this;
    
  if(other.my_obj) // If the other has a non-void object:
  {
    (other.*(other.my_copy_to_other_fnptr))(*this);
  }
  else // If the other has a void object (is INVALID):
  {
    if(my_obj) // If we're currently storing something, delete it and go to an invalid state:
    {
      (this->*my_destroy_me_fnptr)();
    }
    else // If we're storing nothing, go to an invalid state:
    {
      reset();
    }
  }
      
  return *this; 
}
    
//============================================
//  DESTRUCTOR
//============================================
inline UntypedObj::~UntypedObj() 
{ 
  if(my_obj)
    (this->*my_destroy_me_fnptr)(); 
}  

//============================================
//  reset()
//============================================
inline void UntypedObj::reset()
{
  my_obj                 = NULL;
  my_copy_to_other_fnptr = NULL;
  my_destroy_me_fnptr    = NULL;
  my_get_info_fnptr      = NULL;
}

//============================================
//  assignFunctionPtrs()
//============================================
template<typename T> 
inline void UntypedObj::assignFunctionPtrs()
{
  my_copy_to_other_fnptr = copyMeToOther<T>;
  my_destroy_me_fnptr    = destroyMe<T>;
  my_get_info_fnptr      = getInfo<T>;
}
    
//============================================
//  copyMeToOther()
//============================================
template<typename T> 
inline void UntypedObj::copyMeToOther(UntypedObj & other) const
{
  // If the other has a valid object:
  if(other.my_obj)
  {
    // If it's the same type, just assign it:
    if(getInfo() == other.getInfo())
    {
      *((T*)other.my_obj) = *((T*)my_obj);
      return;
    }
    // If it's a different type, delete the old object:
    else
    {
      (other.*(other.my_destroy_me_fnptr))();
    }
  }

  // Now copy over everything:
  other.my_obj                 = new T(*((const T*)(my_obj)));
  other.my_copy_to_other_fnptr = my_copy_to_other_fnptr;   
  other.my_destroy_me_fnptr    = my_destroy_me_fnptr;
  other.my_get_info_fnptr      = my_get_info_fnptr;
}
    
//============================================
//  destroyMe()
//============================================
template<typename T> 
inline void UntypedObj::destroyMe()
{
  if(my_obj)
  {
    delete (T*)my_obj;

    reset();
  }
}

//============================================
//  getInfo()
//============================================
template<typename T> 
inline UTIL::TypeInfo UntypedObj::getInfo() const
{
  return UTIL::TypeInfo::create<T> ();
}

} // End namespace SAGE

#endif
