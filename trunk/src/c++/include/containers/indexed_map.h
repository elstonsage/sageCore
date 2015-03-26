#ifndef SEQUENCE_MAP
#define SEQUENCE_MAP

#include <memory>
#include <algorithm>
#include <iterator>
#include <functional>
#include <stdexcept>
#include <limits>
#include <vector>
#include <map>
#include <cassert>

namespace SAGE 
{

template <class Key, class T, class Compare = std::less<Key>, class Allocator=std::allocator<std::pair<Key,T> > >
class indexed_map
{
public:
    typedef Key                                 key_type;
    typedef std::pair<const Key, T>             value_type;
    typedef Compare                             key_compare;
    typedef Allocator                           allocator_type;
    typedef typename Allocator::reference       reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::pointer         pointer;
    typedef typename Allocator::const_pointer   const_pointer;

private:

    typedef          std::pair<Key,T>                          __element_type;
    typedef          std::vector<__element_type,Allocator>     __element_vector;
    typedef typename __element_vector::iterator                __element_iterator;
    typedef typename __element_vector::const_iterator          __element_const_iterator;

// NOTE: The commented line below was changed to the line immediately below it. I simply
// took out the 'Allocator' part and let the stl use the default allocator for the map.
// Why, you might ask, did I do this? Well, on MacOS at least, the compiler produces tons of
// compiler errors associated with this map. I *think* it's because the Allocator explicitly
// given is for allocating pairs of Key's & T's, NOT pairs of Key's and size_type's.
// --sgross 28 mar 06
//
//    typedef          std::map<Key,size_type,Compare,Allocator> __index_map;
    typedef          std::map<Key,size_type,Compare> __index_map;

    typedef typename __index_map::iterator                     __index_iterator;
    typedef typename __index_map::const_iterator               __index_const_iterator;
    typedef typename __index_map::value_type	               __index_value_type;
    typedef typename __index_map::value_compare                __index_value_compare;

    __index_map      my_index;
    __element_vector my_elements;
    
public:

    class value_compare : public std::binary_function<value_type,value_type,bool> {
    private:
        friend class indexed_map<Key, T, Compare, Allocator>;
        Compare comp;
        value_compare( Compare c ) : comp(c) {}
    public:
        bool operator()(const value_type& x, const value_type& y) {
            return comp(x.first, y.first);
        }
    };

    class indexed_iterator;
    class indexed_const_iterator;
//    typedef std::reverse_iterator<indexed_const_iterator>  indexed_const_reverse_iterator;
//    typedef std::reverse_iterator<indexed_iterator>        indexed_reverse_iterator;

    class ordered_iterator;
    class ordered_const_iterator;
//    typedef std::reverse_iterator<ordered_const_iterator>     ordered_const_reverse_iterator;
//    typedef std::reverse_iterator<ordered_iterator>           ordered_reverse_iterator;

    // - Iterates in index order.
    //
    class indexed_iterator :
                  public std::iterator<std::random_access_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    {
      friend class indexed_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class indexed_const_iterator;

    protected:
      __element_iterator eit;
      indexed_iterator(const __element_iterator &i) : eit(i) {}
      indexed_iterator(const ordered_iterator& i) : eit(i.eit->second) {}
    public:
      indexed_iterator() {}

      bool operator== (const indexed_iterator& x)       const { return eit == x.eit; }
      bool operator== (const indexed_const_iterator& x) const { return eit == x.eit; }
      bool operator!= (const indexed_iterator& x)       const { return eit != x.eit; }
      bool operator!= (const indexed_const_iterator& x) const { return eit != x.eit; }
      bool operator== (const ordered_iterator& x)        const { return *this == x.index(); }
      bool operator== (const ordered_const_iterator& x)  const { return *this == x.index(); }
      bool operator!= (const ordered_iterator& x)        const { return *this != x.index(); }
      bool operator!= (const ordered_const_iterator& x)  const { return *this != x.index(); }

      typedef typename indexed_map<Key,T,Compare,Allocator>::value_type value_type;

      value_type& operator* ()  const { return reinterpret_cast<value_type&>(*eit);   }
      value_type* operator->()  const { return &(operator*()); }
      indexed_iterator& operator++ ()    { ++eit; return *this; }
      indexed_iterator  operator++ (int) { indexed_iterator tmp = *this; ++*this; return tmp; }
      indexed_iterator& operator-- ()    { --eit; return *this; }
      indexed_iterator  operator-- (int) { indexed_iterator tmp = *this; --*this; return tmp; }      

      indexed_iterator& operator+= (int s) { eit+=s; return *this; }
      indexed_iterator  operator+  (int s) { indexed_iterator tmp = *this; tmp+=s; return tmp; }
      indexed_iterator& operator-= (int s) { eit+=s; return *this; }
      indexed_iterator  operator-  (int s) { indexed_iterator tmp = *this; this+=s; return tmp; }      
    };

    class indexed_const_iterator :
                  public std::iterator<std::random_access_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class indexed_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class indexed_iterator;
    protected:
      __element_const_iterator eit;
      indexed_const_iterator(const __element_const_iterator &i)
        : eit(i) {}
    public:
      indexed_const_iterator() {}
      indexed_const_iterator(const indexed_iterator& i) : eit(i) {}
      indexed_const_iterator(const ordered_iterator& i) : eit(i.eit->second) {}
      indexed_const_iterator(const ordered_const_iterator& i) : eit(i.eit->second) {}

      bool operator== (const indexed_iterator& x)       const { return eit == x.eit; }
      bool operator== (const indexed_const_iterator& x) const { return eit == x.eit; }
      bool operator!= (const indexed_iterator& x)       const { return eit != x.eit; }
      bool operator!= (const indexed_const_iterator& x) const { return eit != x.eit; }
      bool operator== (const ordered_iterator& x)          const { return *this == x.index(); }
      bool operator== (const ordered_const_iterator& x)    const { return *this == x.index(); }
      bool operator!= (const ordered_iterator& x)          const { return *this != x.index(); }
      bool operator!= (const ordered_const_iterator& x)    const { return *this != x.index(); }

      typedef typename indexed_map<Key,T,Compare,Allocator>::value_type value_type;

      const value_type& operator* () const { return reinterpret_cast<const value_type&>(*eit); }
      const value_type* operator->() const { return &(operator*());  }
      indexed_const_iterator& operator++ ()    { ++eit; return *this; }
      indexed_const_iterator  operator++ (int) { indexed_const_iterator tmp = *this; ++*this; return tmp; }
      indexed_const_iterator& operator-- ()    { --eit; return *this; }
      indexed_const_iterator  operator-- (int) { indexed_const_iterator tmp = *this; --*this; return tmp; }      

      indexed_const_iterator& operator+= (int s) { eit+=s; return *this; }
      indexed_const_iterator  operator+  (int s) { indexed_iterator tmp = *this; tmp+=s; return tmp; }
      indexed_const_iterator& operator-= (int s) { eit+=s; return *this; }
      indexed_const_iterator  operator-  (int s) { indexed_iterator tmp = *this; this+=s; return tmp; }      
    };

    // - Iterates in MAP order.
    //
    class ordered_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    {
      friend class indexed_map<Key, T, Compare, Allocator>;
      friend class ordered_const_iterator;
      friend class indexed_iterator;
      friend class indexed_const_iterator;
    protected:
      __index_iterator  eit;
      __element_vector* ev;
      ordered_iterator(__element_vector* v, const __index_iterator &i) : eit(i), ev(v) {}

      indexed_iterator       index() const { return indexed_iterator(ev->begin() + eit->second); }

    public:
      ordered_iterator() {}

      bool operator== (const indexed_iterator& x)       const { return index() == x; }
      bool operator== (const indexed_const_iterator& x) const { return index() == x; }
      bool operator!= (const indexed_iterator& x)       const { return index() != x; }
      bool operator!= (const indexed_const_iterator& x) const { return index() != x; }
      bool operator== (const ordered_iterator& x)          const { return eit == x.eit; }
      bool operator== (const ordered_const_iterator& x)    const { return eit == x.eit; }
      bool operator!= (const ordered_iterator& x)          const { return eit != x.eit; }
      bool operator!= (const ordered_const_iterator& x)    const { return eit != x.eit; }

      typedef typename indexed_map<Key,T,Compare,Allocator>::value_type value_type;

      value_type& operator* () const { return reinterpret_cast<value_type&>(*index()); }
      value_type* operator->() const { return &(operator*());  }
      ordered_iterator& operator++ ()   { ++eit; return *this; }
      ordered_iterator  operator++ (int) { ordered_iterator tmp = *this; ++*this; return tmp; }
      ordered_iterator& operator-- ()   { --eit; return *this; }
      ordered_iterator  operator-- (int) { ordered_iterator tmp = *this; --*this; return tmp; }      
    };

    class ordered_const_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    {
      friend class indexed_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class indexed_iterator;
      friend class indexed_const_iterator;
    protected:
      __index_const_iterator  eit;
      const __element_vector* ev;

      ordered_const_iterator(const __element_vector* v, const __index_const_iterator &i) : eit(i), ev(v) {}

      indexed_const_iterator index() const { return indexed_const_iterator(ev->begin() + eit->second); }

    public:
      ordered_const_iterator() {}
      ordered_const_iterator(const ordered_iterator& i) : eit(i) {}

      bool operator== (const indexed_iterator& x)       const { return index() == x; }
      bool operator== (const indexed_const_iterator& x) const { return index() == x; }
      bool operator!= (const indexed_iterator& x)       const { return index() != x; }
      bool operator!= (const indexed_const_iterator& x) const { return index() != x; }
      bool operator== (const ordered_iterator& x)          const { return eit == x.eit; }
      bool operator== (const ordered_const_iterator& x)    const { return eit == x.eit; }
      bool operator!= (const ordered_iterator& x)          const { return eit != x.eit; }
      bool operator!= (const ordered_const_iterator& x)    const { return eit != x.eit; }

      typedef typename indexed_map<Key,T,Compare,Allocator>::value_type value_type;

      const value_type& operator* () const { return reinterpret_cast<const value_type&>(*index()); }
      const value_type* operator->() const { return &(operator*()); }
      ordered_const_iterator& operator++ ()   { ++eit; return *this; }
      ordered_const_iterator operator++ (int) { ordered_const_iterator tmp = *this; ++*this; return tmp; }
      ordered_const_iterator& operator-- ()   { --eit; return *this; }
      ordered_const_iterator operator-- (int) { ordered_const_iterator tmp = *this; --*this; return tmp; }      
    };

    explicit indexed_map(const Compare& comp, const Allocator& allocator) 
               : my_index(comp), my_elements(allocator) { }
    explicit indexed_map(const Compare& comp)
               : my_index(comp) { }
    explicit indexed_map() { }

    template <class InputIterator>
    indexed_map(InputIterator first, InputIterator last,
            const Compare& comp, const Allocator& allocator )
      : my_index(comp), my_elements(allocator)
    {
	insert( first, last );
    }

    template <class InputIterator>
    indexed_map(InputIterator first, InputIterator last, const Compare& comp )
      : my_index(comp), my_elements(Allocator())
    {
	insert( first, last );
    }

    template <class InputIterator>
    indexed_map(InputIterator first, InputIterator last)
      : my_index(Compare()), my_elements(Allocator())
    {
	insert( first, last );
    }

    indexed_map(const indexed_map& x) :
	my_index(x.my_index), my_elements(x.my_elements)
    {}

    allocator_type get_allocator() const {return my_index.get_allocator();}

    indexed_map& operator=(const indexed_map& x)
    {
      my_index = x.my_index;
      my_elements = x.my_elements;
      return *this;
    }

    // ordered iterators:
    ordered_iterator                ordered_begin()         { return ordered_iterator(&my_elements, my_index.begin()); }
    ordered_const_iterator          ordered_begin()   const { return ordered_const_iterator(&my_elements, my_index.begin()); }
    ordered_iterator                ordered_end()           { return ordered_iterator(&my_elements, my_index.end()); }
    ordered_const_iterator          ordered_end()     const { return ordered_const_iterator(&my_elements, my_index.end()); }
//    ordered_reverse_iterator        ordered_rbegin()        { return ordered_reverse_iterator(&my_elements, ordered_end()); }
//    ordered_const_reverse_iterator  ordered_rbegin()  const { return ordered_const_reverse_iterator(&my_elements, ordered_end()); }
//    ordered_reverse_iterator        ordered_rend()          { return ordered_reverse_iterator(&my_elements, ordered_begin()); }
//    ordered_const_reverse_iterator  ordered_rend()    const { return ordered_const_reverse_iterator(&my_elements, ordered_begin()); }

    // indexed iterators:
    indexed_iterator                indexed_begin()         { return indexed_iterator(my_elements.begin()); }
    indexed_const_iterator          indexed_begin()   const { return indexed_const_iterator(my_elements.begin()); }
    indexed_iterator                indexed_end()           { return indexed_iterator(my_elements.end()); }
    indexed_const_iterator          indexed_end()     const { return indexed_const_iterator(my_elements.end()); }
//    indexed_reverse_iterator        indexed_rbegin()        { return indexed_reverse_iterator(indexed_end()); }
//    indexed_const_reverse_iterator  indexed_rbegin()  const { return indexed_const_reverse_iterator(indexed_end()); }
//    indexed_reverse_iterator        indexed_rend()          { return indexed_reverse_iterator(indexed_begin()); }
//    indexed_const_reverse_iterator  indexed_rend()    const { return indexed_const_reverse_iterator(indexed_begin()); }
 
    // capacity:
    bool          empty()     const { return my_index.empty(); }
    size_type     size()      const { return my_index.size();  }
    size_type     max_size()  const { return std::min( my_index.max_size(),
                                                       my_elements.max_size()); }
    // accessors:
    const key_type& key(size_type s) const
    {
      return my_elements[s].first;
    }
    size_type index(const key_type& s) const
    {
      __index_const_iterator i = my_index.find(s);
      
      if(i != my_index.end()) return i->second;
      return (size_type) -1;
    }

    // modifiers:
    T& operator[](const key_type& x)
    {
      value_type y(x,T());
      std::pair<ordered_iterator, bool> ip = insert(y);
      return ip.first->second;
    }

    T& operator[](size_type s)
    {
      return my_elements[s].second;
    }
    const T& operator[](size_type s) const
    {
      return my_elements[s].second;
    }

    void push_back (const value_type& x)
    {
      insert(x);
    }

    void pop_back ()
    {
      __element_iterator i = my_elements.end();
      erase( --i );
    }
    
    std::pair<ordered_iterator, bool> insert(const value_type& x)
    {
      std::pair<__index_iterator, bool> ip = my_index.insert( 
                         __index_value_type(x.first, (size_type) -1 ) );

      if(ip.second)                                      // If we inserted a new
      {                                                  // element, then add it to
        my_elements.push_back(x);                        // the indexed and then
        ip.first->second = my_elements.size() - 1;       // patch it in to the index.
      }
      else                                               // otherwise, just update
      {                                                  // the value in place
        my_elements[ip.first->second].second = x.second; // using the index.
      }

      assert( my_elements.size() == my_index.size() );
      return std::pair<ordered_iterator, bool>(ordered_iterator(&my_elements,ip.first), ip.second);
    }

    ordered_iterator insert(ordered_iterator position, const value_type& x)
    {
      size_type s = size();
      __index_iterator i = my_index.insert( position.eit, 
                          __index_value_type(x.first, (size_type) -1) );

      if( size() != s )                                   // If we inserted a new
      {                                                   // element, then add it to
        my_elements.push_back(x);                         // the indexed and then
        i->second = my_elements.size() - 1;               // patch it in to the index.
      }
      else                                                // otherwise, just update
      {                                                   // the value in place
        my_elements[i->second].second = x.second;         // using the index.
      }

      assert( my_elements.size() == my_index.size() );
      return ordered_iterator(&my_elements,i);
    }
      
    template <class InputIterator> 
    void insert(InputIterator first, InputIterator last)
    {
      if( first == last )
        return;

      ordered_iterator hint = insert( *first++ ).first;
      for( ; first != last; ++first )
        hint = insert( hint, *first );
    }

    void erase(const key_type &x)
    {
      erase(ordered_find(x));
    }
    
    void erase(ordered_iterator position)
    {
      if(position == ordered_end())
        return;
        
      assert( my_elements.size() == my_index.size() );

      const __element_type& last_elt = my_elements[size()-1];
      my_elements[position.eit->second] = last_elt;
      my_index[last_elt.first] = position.eit->second;

      my_elements.resize(size()-1);
      my_index.erase(position.eit);
      assert( my_elements.size() == my_index.size() );
    }

    void erase(indexed_iterator position)
    {
      if(position == indexed_end())
        return;
        
      erase(ordered_find(position->first));
    }

    // Both of the following erases are broken.
/*    
    void erase(ordered_iterator first, ordered_iterator last)
    {
      if( first == last )
        return;

      for( ; first != last; ++first )
        erase( first );
    }

    void erase(indexed_iterator first, indexed_iterator last)
    {
      if( first == last )
        return;

      for( ; first != last; ++first )
        erase( first->first );
    }
*/
    void swap( indexed_map& x )
    {
      // This looks like it's copying the x, not swapping. ??? GCW 001206
      __index_map        temp_map  = x.my_index;
      __element_vector temp_vector = x.my_elements;

      my_index.swap(temp_map);
      my_elements.swap(temp_vector);

      return *this;
    } 

    void clear()
    {
      my_index.clear();
      my_elements.clear();
    }
    
    // _lib.set.observers_ observers:
    key_compare   key_comp() const    { return my_elements.key_comp(); }
    value_compare value_comp() const  { return my_index.value_comp(); }

    // [lib.map.ops] set operations:
    indexed_iterator indexed_find(const key_type& x)
    {
      __index_iterator ei = my_index.find(x);
      if(ei == my_index.end())
        return indexed_end();
      return my_elements.begin() + ei->second;
    }

    indexed_const_iterator indexed_find(const key_type& x) const
    {
      __index_const_iterator ei = my_index.find(x);
      if(ei == my_index.end())
        return indexed_end();
      return my_elements.begin() + ei->second;
    }

    ordered_iterator ordered_find(const key_type& x)
    {
      return ordered_iterator( &my_elements, my_index.find(x) );
    }
    
    ordered_const_iterator ordered_find(const key_type& x) const
    {
      return ordered_const_iterator( &my_elements, my_index.find(x) );
    }

    size_type count(const key_type& x) const
    {
      return my_index.count(x);
    }
    
    // operations:
    ordered_iterator lower_bound(const key_type& x)
    {
      return ordered_iterator( my_index.lower_bound(x) );
    }
    ordered_const_iterator lower_bound(const key_type& x) const
    {
      return ordered_const_iterator( my_index.lower_bound(x) );
    }
    ordered_iterator upper_bound(const key_type& x)
    {
      return ordered_iterator( my_index.upper_bound(x) );
    }
    ordered_const_iterator upper_bound(const key_type& x) const
    {
      return ordered_const_iterator( my_index.upper_bound(x) );
    }
    std::pair<ordered_iterator,ordered_iterator> equal_range(const key_type& x)
    {
      std::pair<__index_iterator, __index_iterator> ip = my_index.equal_range(x);
      return std::pair<ordered_iterator,ordered_iterator>(ip->first,ip->second);
    }

    std::pair<ordered_const_iterator,ordered_const_iterator> equal_range(const key_type& x) const 
    {
      std::pair<__index_const_iterator, __index_const_iterator> ip = my_index.equal_range(x);
      return std::pair<ordered_const_iterator,ordered_const_iterator>(ip->first,ip->second);
    }

};

}

#endif
