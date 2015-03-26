#ifndef CACHE_MAP_H
#define CACHE_MAP_H

#include <functional>
#include <vector>
#include "containers/hash_fun.h"

namespace SAGE {

template <class Key, class T,
          class HashFcn   = hash<Key>,
          class EqualKey  = std::equal_to<Key>,
          class Allocator = std::allocator<T> >
class cache_map
{
private:

  template <class V>
  struct hash_ptr
  {
             hash_ptr() : my_ptr(NULL) { }
    explicit hash_ptr(V* p) : my_ptr(p) { }

    ~hash_ptr() { release(); }

    V*       operator->()       { return  my_ptr; }
    const V* operator->() const { return  my_ptr; }
    V&       operator*()        { return *my_ptr; }
    const V& operator*()  const { return *my_ptr; }

    hash_ptr& operator=(V* v)   { release(v); return *this; }

    template <class Y>
    hash_ptr& operator=(hash_ptr<Y> v)
    {
      release( v.reset() );
      return *this;
    }

    bool     is_set()     const { return my_ptr != NULL; }
    bool     is_null()    const { return my_ptr == NULL; }

    V*       reset(V* x = NULL) 
    { 
      V* t = my_ptr; 
      my_ptr = x; 
      return t; 
    }

    void     release(V* x = NULL)
    {
      if(my_ptr)
        delete(my_ptr);
      my_ptr = x;
    }

  private:
    V*          my_ptr;
  };


public:
  typedef Key                                      key_type;
  typedef T                                        data_type;
  typedef std::pair<const Key, T>                  value_type;
  typedef HashFcn                                  hash_type;
  typedef EqualKey                                 equal_key_type;
  typedef Allocator                                allocator_type;
  typedef typename Allocator::reference            reference;
  typedef typename Allocator::const_reference      const_reference;
  typedef typename Allocator::size_type            size_type;
  typedef typename Allocator::difference_type      difference_type;
  typedef typename Allocator::pointer              pointer;
  typedef typename Allocator::const_pointer        const_pointer;

  hash_type      get_hash_function() const { return my_hash;      }
  equal_key_type get_equal_key()     const { return my_key_eq;    }
  allocator_type get_allocator()     const { return my_allocator; }

  typedef  std::vector< hash_ptr<value_type> >     storage_type;

  typedef  typename storage_type::iterator         iterator;
  typedef  typename storage_type::const_iterator   const_iterator;

public:
  cache_map() : my_storage(769) {}

  explicit cache_map(size_type n) : my_storage(n) { }

  cache_map(size_type n, const hash_type& hf)
    : my_storage(n), my_hash(hf) {}

  cache_map(size_type n, const hash_type& hf, const equal_key_type& eql,
           const allocator_type& a = allocator_type())
    : my_storage(n,a), my_hash(hf), my_key_eq(eql) {}

  template <class InputIterator>
  cache_map(InputIterator f, InputIterator l, size_type n)
    : my_storage(n)
    { insert(f, l); }

  template <class InputIterator>
  cache_map(InputIterator f, InputIterator l, size_type n,
           const hash_type& hf)
    : my_storage(n), my_hash(hf)
    { insert(f, l); }

  template <class InputIterator>
  cache_map(InputIterator f, InputIterator l, size_type n,
           const hash_type& hf, const equal_key_type& eql,
           const allocator_type& a = allocator_type())
    : my_storage(n), my_hash(hf), my_allocator(a)
    { insert(f, l); }

  size_type     size() const { return my_storage.size(); }
  size_type max_size() const { return my_storage.max_size(); }
  bool         empty() const { return my_storage.empty(); }

  void swap(cache_map& cm)
  {
    std::swap(my_storage,   cm.my_storage);
    std::swap(my_allocator, cm.my_allocator);
    std::swap(my_hash,      cm.my_hash);
  }

  template <class K1, class T1, class HF, class EqK, class Al>
  friend bool operator==(const cache_map<K1, T1, HF, EqK, Al>&,
                         const cache_map<K1, T1, HF, EqK, Al>&);

  iterator         begin()           { return my_storage.begin(); }
  iterator         end()             { return my_storage.end(); }
  const_iterator   begin()  const    { return my_storage.begin(); }
  const_iterator   end()    const    { return my_storage.end(); }

  size_t hash(const key_type& k)
  {
    return my_hash(k) % my_storage.size();
  }

  std::pair<iterator,bool> insert(const value_type& obj)
  {
    if( !size() )
      resize(769);

    size_t h = hash(obj.first);

    if( my_storage[h].is_null() )
      my_storage[h]  = new value_type(obj);
    else
    {
      (*my_storage[h]).~value_type();
      new (&*my_storage[h]) value_type(obj);
    }

    return begin() + h;
  }

  iterator find(const key_type& x)
  {
    size_t h = hash(x);
    if( my_storage[h].is_set() && my_key_eq(my_storage[h]->first, x) )
      return begin() + h;
    return end();
  }

  const_iterator find(const key_type& x) const
  {
    return find(x);
  }

  T& operator[](const key_type& x)
  {
    size_t h = hash(x);

    if( my_storage[h].is_null() )
      my_storage[h] = new value_type(x, T());
    else if( !my_key_eq(my_storage[h]->first, x) )
    {
      (*my_storage[h]).~value_type();
      new (&*my_storage[h]) value_type(x, T());
    }
    return my_storage[h]->second;
  }

  template <class InputIterator>
  void insert(InputIterator f, InputIterator l)
  {
    insert(f,l);
  }

  size_type erase(const key_type& key)
  {
    size_t h = hash(key);
    if(my_storage[h].is_set())
      my_storage[h].release();
  }

  void erase(iterator it)
  {
    it = value_type();
  }

  void erase(iterator f, iterator l)
  {
    erase(f, l);
  }

  void clear()
  {
    erase(begin(), end());
  }

  void resize(size_type n)
  {
    if(n == size())
      return;

    storage_type s(n);
    for(size_t i = 0; i < my_storage.size(); ++i)
      if( my_storage[i].is_set() )
        s[ hash(my_storage[i]->first) ] = my_storage[i];

    std::swap(s, my_storage);
  }

private:

  storage_type   my_storage;
  allocator_type my_allocator;
  hash_type      my_hash;
  equal_key_type    my_key_eq;
};

template <class Key, class T, class HashFcn, class EqlKey, class Allocator>
inline void
swap(cache_map<Key,T,HashFcn,EqlKey,Allocator>& hm1,
     cache_map<Key,T,HashFcn,EqlKey,Allocator>& hm2)
{
  hm1.swap(hm2);
}

}

#if 0
namespace std
{

template <class Key, class T, class HashFn,  class EqKey, class Allocator>
class insert_iterator<SAGE::cache_map<Key, T, HashFn, EqKey, Allocator> > {
protected:
  typedef SAGE::cache_map<Key, T, HashFn, EqKey, Allocator> Container;
  Container* container;
public:
  typedef Container          container_type;
  typedef std::output_iterator_tag iterator_category;
  typedef void                value_type;
  typedef void                difference_type;
  typedef void                pointer;
  typedef void                reference;

  insert_iterator(Container& x) : container(&x) {}
  insert_iterator(Container& x, typename Container::iterator)
    : container(&x) {}
  insert_iterator<Container>&
  operator=(const typename Container::value_type& value) {
    container->insert(value);
    return *this;
  }
  insert_iterator<Container>& operator*()     { return *this; }
  insert_iterator<Container>& operator++()    { return *this; }
  insert_iterator<Container>& operator++(int) { return *this; }
};

}
#endif

#endif
