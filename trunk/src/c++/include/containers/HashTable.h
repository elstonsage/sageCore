#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <map>

namespace SAGE       {
namespace CONTAINERS {

/**
  * The purpose of a HashTable is to provide a nice-looking object-oriented
  * interface for an STL map. Given a KeyType and a ValueType, the HashTable
  * lets you manage values associated with keys.
  */
template<class KeyType, class ValueType>
class HashTable
{
  public:
  
  /// @name Constructors
  //@{
  
    /**
      * Constructor.
      */
    HashTable() { my_map.clear(); }
    
    /**
      * Copy constructor.
      */
    HashTable(const HashTable & other) { my_map = other.my_map; }
    
  //@}

  /// @name Operators
  //@{
  
    /**
      * Operator=
      */
    HashTable & operator= (const HashTable & other) { if(this == &other) return *this; my_map = other.my_map; return *this; }
  
  //@}

  /// @name Adding / setting values
  //@{

    /**
      * Adds the indicated value with the specified key.
      * \param key The key associated with the value
      * \param value The value to add
      * \retval 0 If successful
      * \retval 1 If not successful (if key already exists, for instance)
      */
    int addValue(const KeyType & key, const ValueType & value);
    
    /**
      * Sets the value associated with the given key.
      * 
      * If the key does not exist, setValue creates an entry for that key.
      * \param key The key associated with the value
      * \param value The value to associate with the given key
      */
    void setValue(const KeyType & key, const ValueType & value);
    
  //@}
    
  /// @name Getting values
  //@{
  
    /**
      * Returns the number of key-value pairs stored in this object.
      */
    int size() const;

    /**
      * Indicates whether or not there is an entry for the given key.
      * \param key The key in question
      * \returns true If the key does exist
      * \returns false If the key does not exist
      */
    bool doesKeyExist(const KeyType & key) const;

    /**
      * Returns the value associated with the given key.
      * \param key The key associated with the requested value
      * \returns The value associated with the given key
      */
    ValueType & getValue(const KeyType & key);

    /**
      * Returns the value associated with the given key.
      * \param key The key associated with the requested value
      * \returns The value associated with the given key
      */
    const ValueType & getValue(const KeyType & key) const;
    
    /**
      * Returns one of the values stored in this object.
      *
      * This function is useful if there is some aspect 
      * of your values that is constant across all values.
      * In that case, fetching any value is sufficient for
      * accessing that constant information.
      */
    const ValueType & getAnyValue() const;

  //@}

  private:
  
    std::map<KeyType, ValueType> my_map;
};

//============================
//  INLINE FUNCTIONS
//============================

//============================
//
//  size
//
//============================
template<class KeyType, class ValueType>
inline int 
HashTable<KeyType, ValueType>::size() const
{
  return my_map.size();
}

//============================
//
//  doesKeyExist
//
//============================
template<class KeyType, class ValueType>
inline bool 
HashTable<KeyType, ValueType>::doesKeyExist(const KeyType & key) const
{
  return my_map.find(key) != my_map.end();
}

//============================
//
//  addValue
//
//============================
template<class KeyType, class ValueType>
inline int 
HashTable<KeyType, ValueType>::addValue(const KeyType & key, const ValueType & value)
{
  if(doesKeyExist(key) == 0)
  {
    my_map[key] = value;
    return 0;
  }
  else
  {
    return 1;
  }
}

//============================
//
//  setValue
//
//============================
template<class KeyType, class ValueType>
inline void 
HashTable<KeyType, ValueType>::setValue(const KeyType & key, const ValueType & value)
{
  my_map[key] = value;
}

//============================
//
//  getValue NON-CONST
//
//============================
template<class KeyType, class ValueType>
inline ValueType & 
HashTable<KeyType, ValueType>::getValue(const KeyType & key)
{
  assert(doesKeyExist(key));
  return my_map[key];
}

//============================
//
//  getValue CONST
//
//============================
template<class KeyType, class ValueType>
inline const ValueType & 
HashTable<KeyType, ValueType>::getValue(const KeyType & key) const
{
  assert(doesKeyExist(key));
  return my_map.find(key)->second;
}

//============================
//
//  getAnyValue
//
//============================
template<class KeyType, class ValueType>
inline const ValueType & 
HashTable<KeyType, ValueType>::getAnyValue() const
{
  assert(my_map.size());
  return my_map.begin()->second;
}

} // End namespace CONTAINERS
} // End namespace SAGE

#endif
