#ifndef MSIPL_VECTOR_H
#define MSIPL_VECTOR_H

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

template <class Key, class Compare = std::less<Key>, class Allocator=std::allocator<Key> >
class sequence_set
{
    typedef std::list<const Key*>                   __sequence_list;
    typedef __sequence_list::iterator               __sequence_iterator;
    typedef __sequence_list::const_iterator         __sequence_const_iterator;
    typedef __sequence_list::reverse_iterator       __sequence_reverse_iterator;
    typedef __sequence_list::const_reverse_iterator __sequence_const_reverse_iterator;

    typedef std::map<Key,__sequence_iterator,Compare,Allocator> __element_map;
    typedef __element_map::iterator                 __element_iterator;
    typedef __element_map::const_iterator           __element_const_iterator;
    typedef __element_map::reverse_iterator         __element_reverse_iterator;
    typedef __element_map::const_reverse_iterator   __element_const_reverse_iterator;
    typedef __element_map::value_type	            __element_value_type;
    typedef __element_map::value_compare            __element_value_compare;

    __sequence_list my_sequence;
    __element_map   my_elements;
    
public:

    typedef Key                                 key_type;
    typedef Key                                 value_type;
    typedef Compare                             key_compare;
    typedef Compare                             value_compare;
    typedef Allocator                           allocator_type;
    typedef typename Allocator::reference       reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::size_type       size_type;
    typedef typename Allocator::difference_type difference_type;
    typedef typename Allocator::pointer         pointer;
    typedef typename Allocator::const_pointer   const_pointer;

    class sequence_iterator;
    class sequence_const_iterator;

    typedef std::reverse_iterator<sequence_const_iterator>  sequence_const_reverse_iterator;
    typedef std::reverse_iterator<sequence_iterator>        sequence_reverse_iterator;

    class ordered_iterator;
    class ordered_const_iterator;
    typedef std::reverse_iterator<ordered_const_iterator>     ordered_const_reverse_iterator;
    typedef std::reverse_iterator<ordered_iterator>           ordered_reverse_iterator;

    class sequence_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_set<Key, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class sequence_const_iterator;

    protected:
      __sequence_iterator eit;
      sequence_iterator(const __sequence_iterator &i) : eit(i) {}
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

      const_reference operator* ()   const { return **eit;   }
      const_pointer   operator-> ()  const { return &(operator*()); }
      sequence_iterator& operator++ ()    { ++eit; return *this; }
      sequence_iterator  operator++ (int) { sequence_iterator tmp = *this; ++*this; return tmp; }
      sequence_iterator& operator-- ()    { --eit; return *this; }
      sequence_iterator  operator-- (int) { sequence_iterator tmp = *this; --*this; return tmp; }      
    };

    class sequence_const_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_set<Key, Compare, Allocator>;
      friend class ordered_iterator;
      friend class ordered_const_iterator;
      friend class sequence_iterator;
    protected:
      __sequence_const_iterator eit;
      sequence_const_iterator(const __sequence_const_iterator &i) : eit(i) {}
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

      const_reference operator* () const { return **eit; }
      const_pointer operator-> ()  const { return &(operator*());  }
      sequence_const_iterator& operator++ ()    { ++eit; return *this; }
      sequence_const_iterator  operator++ (int) { sequence_const_iterator tmp = *this; ++*this; return tmp; }
      sequence_const_iterator& operator-- ()    { --eit; return *this; }
      sequence_const_iterator  operator-- (int) { sequence_const_iterator tmp = *this; --*this; return tmp; }      
    };

    class ordered_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_set<Key, Compare, Allocator>;
      friend class ordered_const_iterator;
      friend class sequence_iterator;
      friend class sequence_const_iterator;
    protected:
      __element_iterator eit;
      ordered_iterator(const __element_iterator &i) : eit(i) {}
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

      const_reference operator* () const { return eit->first; }
      const_pointer operator-> ()  const { return &(operator*());  }
      ordered_iterator& operator++ ()    { ++eit; return *this; }
      ordered_iterator  operator++ (int) { ordered_iterator tmp = *this; ++*this; return tmp; }
      ordered_iterator& operator-- ()    { --eit; return *this; }
      ordered_iterator  operator-- (int) { ordered_iterator tmp = *this; --*this; return tmp; }      
    };

    class ordered_const_iterator :
                  public std::iterator<std::bidirectional_iterator_tag, key_type,
                                       difference_type, pointer, reference>
    {
      friend class sequence_set<Key, Compare, Allocator>;
      friend class ordered_iterator;
      friend class sequence_iterator;
      friend class sequence_const_iterator;
    protected:
      __element_iterator eit;
      ordered_const_iterator(const __element_const_iterator &i) : eit(i) {}
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

      const_reference operator* () const { return eit->first; }
      const_pointer operator-> ()  const { return &(operator*()); }
      ordered_const_iterator& operator++ ()   { ++eit; return *this; }
      ordered_const_iterator operator++ (int) { ordered_const_iterator tmp = *this; ++*this; return tmp; }
      ordered_const_iterator& operator-- ()   { --eit; return *this; }
      ordered_const_iterator operator-- (int) { ordered_const_iterator tmp = *this; --*this; return tmp; }      
    };

    explicit sequence_set(const Compare& comp, const Allocator& allocator) 
               : my_elements(comp, allocator) { }
    explicit sequence_set(const Compare& comp)
               : my_elements(comp) { }
    explicit sequence_set() { }

    template <class InputIterator>
    sequence_set(InputIterator first, InputIterator last,
            const Compare& comp, const Allocator& allocator )
      : my_elements(comp,allocator)
    {
	insert( first, last );
    }

    template <class InputIterator>
    sequence_set(InputIterator first, InputIterator last, const Compare& comp )
      : my_elements(comp)
    {
	insert( first, last );
    }

    template <class InputIterator>
    sequence_set(InputIterator first, InputIterator last)
    {
	insert( first, last );
    }

    sequence_set(const sequence_set& x) :
	my_elements(x.my_elements), my_sequence(x.my_sequence)
    {}

    allocator_type get_allocator() const {return value_allocator;}

    sequence_set& operator=(const sequence_set& x)
    {
      my_elements = x.my_elements;
      my_sequence = x.my_sequence;
      return *this;
    }

    // ordered iterators:
    ordered_iterator                ordered_begin()         { return ordered_iterator(my_elements.begin()); }
    ordered_const_iterator          ordered_begin()   const { return ordered_const_iterator(my_elements.begin()); }
    ordered_iterator                ordered_end()           { return ordered_iterator(my_elements.end()); }
    ordered_const_iterator          ordered_end()     const { return ordered_const_iterator(my_elements.end()); }
    ordered_reverse_iterator        ordered_rbegin()        { return ordered_reverse_iterator(ordered_end()); }
    ordered_const_reverse_iterator  ordered_rbegin()  const { return ordered_const_reverse_iterator(ordered_end()); }
    ordered_reverse_iterator        ordered_rend()          { return ordered_reverse_iterator(ordered_begin()); }
    ordered_const_reverse_iterator  ordered_rend()    const { return ordered_const_reverse_iterator(ordered_begin()); }

    // sequence iterators:
    sequence_iterator                sequence_begin()         { return sequence_iterator(my_sequence.begin()); }
    sequence_const_iterator          sequence_begin()   const { return sequence_const_iterator(my_sequence.begin()); }
    sequence_iterator                sequence_end()           { return sequence_iterator(my_sequence.end()); }
    sequence_const_iterator          sequence_end()     const { return sequence_const_iterator(my_sequence.end()); }
    sequence_reverse_iterator        sequence_rbegin()        { return sequence_reverse_iterator(sequence_end()); }
    sequence_const_reverse_iterator  sequence_rbegin()  const { return sequence_const_reverse_iterator(sequence_end()); }
    sequence_reverse_iterator        sequence_rend()          { return sequence_reverse_iterator(sequence_begin()); }
    sequence_const_reverse_iterator  sequence_rend()    const { return sequence_const_reverse_iterator(sequence_begin()); }
 
    // capacity:
    bool          empty()     const { return my_elements.empty(); }
    size_type     size()      const { return my_elements.size();  }
    size_type     max_size()  const { return std::min( my_elements.max_size(),
                                                       my_sequence.max_size()); }
    // modifiers:
    void push_back (const key_type& x)
    {
      insert(x);
    }

    void pop_back ()
    {
      erase( --my_sequence.end() );
    }
    
    std::pair<ordered_iterator, bool> insert(const key_type& x)
    {
      std::pair<__element_iterator, bool> ip = my_elements.insert( 
                                 __element_value_type(x, __sequence_iterator()) );
      if(ip.second)
      {
        my_sequence.push_back( &ip.first->first );
        ip.first->second = --my_sequence.end();
      }
      assert( my_sequence.size() == my_elements.size() );
      return std::pair<ordered_iterator, bool>(ip.first, ip.second);
    }

    ordered_iterator insert(ordered_iterator position, const key_type& x)
    {
      size_type s = size();
      ordered_iterator ip = my_elements.insert( position, 
                                __element_value_type(x, __sequence_iterator()) );

      if( size() != s )
      {
         my_sequence.push_back( &ip.first->first );
         ip.first->second = --my_sequence.end();
      }

      assert( my_sequence.size() == my_elements.size() );
      return ip;
    }
      
    template <class InputIterator> 
    void insert(InputIterator first, InputIterator last)
    {
      if( first == last )
        return;

      for( ; first != last; ++first )
        insert( *first );
    }

    void erase(const key_type &x)
    {
      ordered_iterator oi = my_elements.find(x);
      erase(oi);
    }
    
    void erase(ordered_iterator position)
    {
      if(position == ordered_end())
        return;
        
      assert( my_sequence.size() == my_elements.size() );
      my_sequence.erase(position.eit->second);
      my_elements.erase(position.eit);
      assert( my_sequence.size() == my_elements.size() );
    }

    void erase(sequence_iterator position)
    {
      if(position == sequence_end())
        return;
        
      ordered_iterator oi = my_elements.find(*position);
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

    void swap( sequence_set& x )
    {
      __element_map   temp_map  = x.my_elements;
      __sequence_list temp_list = x.my_sequence;

      my_elements.swap(temp_map);
      my_sequence.swap(temp_list);

      return *this;
    } 

    void clear()
    {
      my_elements.clear();
      my_sequence.clear();
    }
    
    // _lib.set.observers_ observers:
    key_compare   key_comp() const    { return my_elments.key_comp(); }
    value_compare value_comp() const  { return my_elements.value_comp(); }

    // [lib.set.ops] set operations:
    sequence_iterator     sequence_find(const key_type& x)
    {
      __element_iterator ei = my_elements.find(x);
      if(ei == my_elements.end())
        return sequence_end();
      return ei->second;
    }

    sequence_const_iterator sequence_find(const key_type& x) const
    {
      __element_iterator ei = my_elements.find(x);
      if(ei == my_elements.end())
        return sequence_end();
      return ei->second;
    }

    ordered_iterator ordered_find(const key_type& x)
    {
      return ordered_iterator( my_elements.find(x) );
    }
    
    ordered_const_iterator ordered_find(const key_type& x) const
    {
      return ordered_const_iterator( my_elemements.find(x) );
    }

    ordered_iterator find(const key_type& x)
    {
      return ordered_find(x);
    }
    
    ordered_const_iterator find(const key_type& x) const
    {
      return ordered_find(x);
    }

    size_type count(const key_type& x) const
    {
      return my_elements.count(x);
    }
    
    // operations:
    ordered_iterator lower_bound(const key_type& x)
    {
      return ordered_iterator( my_elements.lower_bound(x) );
    }
    ordered_const_iterator lower_bound(const key_type& x) const
    {
      return ordered_const_iterator( my_elements.lower_bound(x) );
    }
    ordered_iterator upper_bound(const key_type& x)
    {
      return ordered_iterator( my_elements.upper_bound(x) );
    }
    ordered_const_iterator upper_bound(const key_type& x) const
    {
      return ordered_const_iterator( my_elements.upper_bound(x) );
    }

    std::pair<ordered_iterator,ordered_iterator> equal_range(const key_type& x)
    {
      std::pair<__element_iterator, __element_iterator> ip = my_elements.equal_range(x);
      return std::pair<ordered_iterator,ordered_iterator>(ip->first,ip->second);
    }

    std::pair<ordered_const_iterator,ordered_const_iterator> equal_range(const key_type& x) const 
    {
      std::pair<__element_const_iterator, __element_const_iterator> ip = my_elements.equal_range(x);
      return std::pair<ordered_const_iterator,ordered_const_iterator>(ip->first,ip->second);
    }

};

}

#endif
