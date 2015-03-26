#ifndef SEQUENCE_MAP
#define SEQUENCE_MAP

#include <memory>
#include <algorithm>
#include <iterator>
#include <functional>
#include <stdexcept>
#include <limits>
#include <list>
#include <map>
#include <cassert>

namespace SAGE 
{

template <class Key, class T, class Compare = std::less<Key>, class Allocator=std::allocator<std::pair<Key,T> > >
class sequence_map
{
    typedef std::pair<Key,T>                       __element_type;
    typedef std::list<__element_type,Allocator>    __element_list;
    typedef typename __element_list::iterator               __element_iterator;
    typedef typename __element_list::const_iterator         __element_const_iterator;
    typedef typename __element_list::reverse_iterator       __element_reverse_iterator;
    typedef typename __element_list::const_reverse_iterator __element_const_reverse_iterator;

// To understand why the line below was changed, see indexed_map.h for a similar example - sgross 29 mar 06
//    typedef std::map<Key,__element_iterator,Compare,Allocator> __index_map;
    typedef std::map<Key,__element_iterator,Compare> __index_map;

    typedef typename __index_map::iterator                 __index_iterator;
    typedef typename __index_map::const_iterator           __index_const_iterator;
    typedef typename __index_map::reverse_iterator         __index_reverse_iterator;
    typedef typename __index_map::const_reverse_iterator   __index_const_reverse_iterator;
    typedef typename __index_map::value_type	          __index_value_type;
    typedef typename __index_map::value_compare            __index_value_compare;

    __element_list my_elements;
    __index_map my_index;
    
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

    class value_compare : public std::binary_function<value_type,value_type,bool> {
    private:
        friend class sequence_map<Key, T, Compare, Allocator>;
        Compare comp;
        value_compare( Compare c ) : comp(c) {}
    public:
        bool operator()(const value_type& x, const value_type& y) {
            return comp(x.first, y.first);
        }
    };

    class sequence_iterator;
    class sequence_const_iterator;
    typedef std::reverse_iterator<sequence_const_iterator>  sequence_const_reverse_iterator;
    typedef std::reverse_iterator<sequence_iterator>        sequence_reverse_iterator;

    class ordered_iterator;
    class ordered_const_iterator;
    typedef std::reverse_iterator<ordered_const_iterator>     ordered_const_reverse_iterator;
    typedef std::reverse_iterator<ordered_iterator>           ordered_reverse_iterator;

    class sequence_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class sequence_const_iterator;

    protected:
      __element_iterator eit;
      sequence_iterator(const __element_iterator &i) : eit(i) {}
      sequence_iterator(const ordered_iterator& i) : eit(i.eit->second) {}
    public:
      sequence_iterator() {}

      bool operator== (const sequence_iterator& x)       const { return eit == x.eit; }
      bool operator== (const sequence_const_iterator& x) const { return eit == x.eit; }
      bool operator!= (const sequence_iterator& x)       const { return eit != x.eit; }
      bool operator!= (const sequence_const_iterator& x) const { return eit != x.eit; }
      bool operator== (const ordered_iterator& x)        const { return eit == x.eit->second; }
      bool operator== (const ordered_const_iterator& x)  const { return eit == x.eit->second; }
      bool operator!= (const ordered_iterator& x)        const { return eit != x.eit->second; }
      bool operator!= (const ordered_const_iterator& x)  const { return eit != x.eit->second; }

      typename sequence_iterator::value_type& operator* ()  const { return reinterpret_cast<typename sequence_iterator::value_type&>(*eit);   }
      typename sequence_iterator::value_type* operator->()  const { return &(operator*()); }
      sequence_iterator& operator++ ()    { ++eit; return *this; }
      sequence_iterator  operator++ (int) { sequence_iterator tmp = *this; ++*this; return tmp; }
      sequence_iterator& operator-- ()    { --eit; return *this; }
      sequence_iterator  operator-- (int) { sequence_iterator tmp = *this; --*this; return tmp; }      
    };

    class sequence_const_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class sequence_iterator;
    protected:
      __element_const_iterator eit;
      sequence_const_iterator(const __element_const_iterator &i) : eit(i) {}
    public:
      sequence_const_iterator() {}
      sequence_const_iterator(const sequence_iterator& i) : eit(i) {}
      sequence_const_iterator(const ordered_iterator& i) : eit(i.eit->second) {}
      sequence_const_iterator(const ordered_const_iterator& i) : eit(i.eit->second) {}

      bool operator== (const sequence_iterator& x)       const { return eit == x.eit; }
      bool operator== (const sequence_const_iterator& x) const { return eit == x.eit; }
      bool operator!= (const sequence_iterator& x)       const { return eit != x.eit; }
      bool operator!= (const sequence_const_iterator& x) const { return eit != x.eit; }
      bool operator== (const ordered_iterator& x)          const { return eit->second == x.eit; }
      bool operator== (const ordered_const_iterator& x)    const { return eit->second == x.eit; }
      bool operator!= (const ordered_iterator& x)          const { return eit->second != x.eit; }
      bool operator!= (const ordered_const_iterator& x)    const { return eit->second != x.eit; }

      const typename sequence_iterator::value_type& operator* () const { return reinterpret_cast<const typename sequence_iterator::value_type&>(*eit); }
      const typename sequence_iterator::value_type* operator->() const { return &(operator*());  }
      sequence_const_iterator& operator++ ()    { ++eit; return *this; }
      sequence_const_iterator  operator++ (int) { sequence_const_iterator tmp = *this; ++*this; return tmp; }
      sequence_const_iterator& operator-- ()    { --eit; return *this; }
      sequence_const_iterator  operator-- (int) { sequence_const_iterator tmp = *this; --*this; return tmp; }      
    };

    class ordered_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_map<Key, T, Compare, Allocator>;
      friend class ordered_const_iterator;
      friend class sequence_iterator;
      friend class sequence_const_iterator;
    protected:
      __index_iterator eit;
      ordered_iterator(const __index_iterator &i) : eit(i) {}
    public:
      ordered_iterator() {}

      bool operator== (const sequence_iterator& x)       const { return eit->second == x.eit; }
      bool operator== (const sequence_const_iterator& x) const { return eit->second == x.eit; }
      bool operator!= (const sequence_iterator& x)       const { return eit->second != x.eit; }
      bool operator!= (const sequence_const_iterator& x) const { return eit->second != x.eit; }
      bool operator== (const ordered_iterator& x)          const { return eit == x.eit; }
      bool operator== (const ordered_const_iterator& x)    const { return eit == x.eit; }
      bool operator!= (const ordered_iterator& x)          const { return eit != x.eit; }
      bool operator!= (const ordered_const_iterator& x)    const { return eit != x.eit; }
      typename sequence_iterator::value_type& operator* () const { return reinterpret_cast<typename sequence_iterator::value_type&>(*eit->second); }
      typename sequence_iterator::value_type* operator->() const { return &(operator*());  }
      ordered_iterator& operator++ ()   { ++eit; return *this; }
      ordered_iterator  operator++ (int) { ordered_iterator tmp = *this; ++*this; return tmp; }
      ordered_iterator& operator-- ()   { --eit; return *this; }
      ordered_iterator  operator-- (int) { ordered_iterator tmp = *this; --*this; return tmp; }      
    };

    class ordered_const_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, value_type,
                                       difference_type, pointer, reference>
    { 	
      friend class sequence_map<Key, T, Compare, Allocator>;
      friend class ordered_iterator;
      friend class sequence_iterator;
      friend class sequence_const_iterator;
    protected:
      __index_iterator eit;
      ordered_const_iterator(const __index_const_iterator &i) : eit(i) {}
    public:
      ordered_const_iterator() {}
      ordered_const_iterator(const ordered_iterator& i) : eit(i) {}

      bool operator== (const sequence_iterator& x)       const { return eit->second == x.eit; }
      bool operator== (const sequence_const_iterator& x) const { return eit->second == x.eit; }
      bool operator!= (const sequence_iterator& x)       const { return eit->second != x.eit; }
      bool operator!= (const sequence_const_iterator& x) const { return eit->second != x.eit; }
      bool operator== (const ordered_iterator& x)          const { return eit == x.eit; }
      bool operator== (const ordered_const_iterator& x)    const { return eit == x.eit; }
      bool operator!= (const ordered_iterator& x)          const { return eit != x.eit; }
      bool operator!= (const ordered_const_iterator& x)    const { return eit != x.eit; }

      const typename sequence_iterator::value_type& operator* () const { return reinterpret_cast<const typename sequence_iterator::value_type&>(*eit->second); }
      const typename sequence_iterator::value_type* operator->() const { return &(operator*()); }
      ordered_const_iterator& operator++ ()   { ++eit; return *this; }
      ordered_const_iterator operator++ (int) { ordered_const_iterator tmp = *this; ++*this; return tmp; }
      ordered_const_iterator& operator-- ()   { --eit; return *this; }
      ordered_const_iterator operator-- (int) { ordered_const_iterator tmp = *this; --*this; return tmp; }      
    };

    explicit sequence_map(const Compare& comp, const Allocator& allocator) 
               : my_index(comp), my_elements(allocator) { }
    explicit sequence_map(const Compare& comp)
               : my_index(comp) { }
    explicit sequence_map() { }

    template <class InputIterator>
    sequence_map(InputIterator first, InputIterator last,
            const Compare& comp, const Allocator& allocator )
      : my_index(comp), my_elements(allocator)
    {
	insert( first, last );
    }

    template <class InputIterator>
    sequence_map(InputIterator first, InputIterator last, const Compare& comp )
      : my_index(comp), my_elements(Allocator())
    {
	insert( first, last );
    }

    template <class InputIterator>
    sequence_map(InputIterator first, InputIterator last)
      : my_index(Compare()), my_elements(Allocator())
    {
	insert( first, last );
    }

    sequence_map(const sequence_map& x) :
	my_index(x.my_index), my_elements(x.my_elements)
    {}

// Temporarily commented out: sgross 1 Nov 2005
//    allocator_type get_allocator() const {return value_allocator;}

    sequence_map& operator=(const sequence_map& x)
    {
      my_index = x.my_index;
      my_elements = x.my_elements;
      return *this;
    }

    // ordered iterators:
    ordered_iterator                ordered_begin()         { return ordered_iterator(my_index.begin()); }
    ordered_const_iterator          ordered_begin()   const { return ordered_const_iterator(my_index.begin()); }
    ordered_iterator                ordered_end()           { return ordered_iterator(my_index.end()); }
    ordered_const_iterator          ordered_end()     const { return ordered_const_iterator(my_index.end()); }
    ordered_reverse_iterator        ordered_rbegin()        { return ordered_reverse_iterator(ordered_end()); }
    ordered_const_reverse_iterator  ordered_rbegin()  const { return ordered_const_reverse_iterator(ordered_end()); }
    ordered_reverse_iterator        ordered_rend()          { return ordered_reverse_iterator(ordered_begin()); }
    ordered_const_reverse_iterator  ordered_rend()    const { return ordered_const_reverse_iterator(ordered_begin()); }

    // sequence iterators:
    sequence_iterator                sequence_begin()         { return sequence_iterator(my_elements.begin()); }
    sequence_const_iterator          sequence_begin()   const { return sequence_const_iterator(my_elements.begin()); }
    sequence_iterator                sequence_end()           { return sequence_iterator(my_elements.end()); }
    sequence_const_iterator          sequence_end()     const { return sequence_const_iterator(my_elements.end()); }
    sequence_reverse_iterator        sequence_rbegin()        { return sequence_reverse_iterator(sequence_end()); }
    sequence_const_reverse_iterator  sequence_rbegin()  const { return sequence_const_reverse_iterator(sequence_end()); }
    sequence_reverse_iterator        sequence_rend()          { return sequence_reverse_iterator(sequence_begin()); }
    sequence_const_reverse_iterator  sequence_rend()    const { return sequence_const_reverse_iterator(sequence_begin()); }
 
    // capacity:
    bool          empty()     const { return my_index.empty(); }
    size_type     size()      const { return my_index.size();  }
    size_type     max_size()  const { return std::min( my_index.max_size(),
                                                       my_elements.max_size()); }
    // modifiers:
    T& operator[](const key_type& x)
    {
      value_type y(x,T());
      std::pair<ordered_iterator, bool> ip = insert(y);
      return ip.first->second;
    }

    void push_back (const value_type& x)
    {
      insert(x);
    }

    void pop_back ()
    {
      erase( --my_elements.end() );
    }
    
    std::pair<ordered_iterator, bool> insert(const value_type& x)
    {
      std::pair<__index_iterator, bool> ip = my_index.insert( 
                         __index_value_type(x.first, __element_iterator()) );

      if(ip.second)                               // If we inserted a new
      {                                           // element, then add it to
        my_elements.push_back(x);                 // the sequence and then
        ip.first->second = --my_elements.end();   // patch it in to the index.
      }
      else                                        // otherwise, just update
      {                                           // the value in place
        ip.first->second->second = x.second;      // using the index.
      }

      assert( my_elements.size() == my_index.size() );
      return std::pair<ordered_iterator, bool>(ip.first, ip.second);
    }

    ordered_iterator insert(ordered_iterator position, const value_type& x)
    {
      size_type s = size();
      ordered_iterator ip = my_index.insert( position, 
                          __index_value_type(x.first, __element_iterator()) );

      if( size() != s )                           // If we inserted a new
      {                                           // element, then add it to
        my_elements.push_back(x);                 // the sequence and then
        ip.first->second = --my_elements.end();   // patch it in to the index.
      }
      else                                        // otherwise, just update
      {                                           // the value in place
        ip.first->second->second = x.second;      // using the index.
      }

      assert( my_elements.size() == my_index.size() );
      return ip;
    }
      
    template <class InputIterator> 
    void insert(InputIterator first, InputIterator last)
    {
      if( first == last )
        return;

      ordered_iterator hint = insert( *first++ );
      for( ; first != last; ++first )
        hint = insert( hint, *first );
    }

    void erase(const key_type &x)
    {
      ordered_iterator oi = my_index.find(x);
      erase(oi);
    }
    
    void erase(ordered_iterator position)
    {
      if(position == ordered_end())
        return;
        
      assert( my_elements.size() == my_index.size() );
      my_elements.erase(position.eit->second);
      my_index.erase(position.eit);
      assert( my_elements.size() == my_index.size() );
    }

    void erase(sequence_iterator position)
    {
      if(position == sequence_end())
        return;
        
      ordered_iterator oi = my_index.find(position->first);
      erase(oi);
    }

    void erase(ordered_iterator first, ordered_iterator last)
    {
      if( first == last )
        return;

      for( ; first != last; ++first )
        erase( *first );
    }

    void erase(sequence_iterator first, sequence_iterator last)
    {
      if( first == last )
        return;

      for( ; first != last; ++first )
        erase( *first );
    }

    void swap( sequence_map& x )
    {
      __index_map    temp_map  = x.my_index;
      __element_list temp_list = x.my_elements;

      my_index.swap(temp_map);
      my_elements.swap(temp_list);

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
    sequence_iterator sequence_find(const key_type& x)
    {
      __index_iterator ei = my_index.find(x);
      if(ei == my_index.end())
        return sequence_end();
      return ei->second;
    }

    sequence_const_iterator sequence_find(const key_type& x) const
    {
      __index_iterator ei = my_index.find(x);
      if(ei == my_index.end())
        return sequence_end();
      return ei->second;
    }

    ordered_iterator ordered_find(const key_type& x)
    {
      return ordered_iterator( my_index.find(x) );
    }
    
    ordered_const_iterator ordered_find(const key_type& x) const
    {
      return ordered_const_iterator( my_elements.find(x) );
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
