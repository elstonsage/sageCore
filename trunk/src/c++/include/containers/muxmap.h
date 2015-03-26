#ifndef MUXMAP_H
#define MUXMAP_H

#include <functional>
#include <iterator>
#include <set>
#include <map>

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     MultiplexMap                                                 ~
// ~                                                                         ~
// ~ Purpose:   Define Multiplex Map container class.                        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <
           class K,
           class T,
           class TypeCompare   = std::less<K>,
           class ValueCompare  = std::less<T>,
           class Container     = std::set<T, ValueCompare>
         >
class MultiplexMap
{
public:
  class iterator;
  class const_iterator;
  friend class iterator;
  friend class const_iterator;

  typedef  MultiplexMap<K,T,TypeCompare,
                        ValueCompare,Container>        self_type;
  
  typedef  typename Container::size_type               size_type;
  typedef  typename Container::difference_type         difference_type;
  typedef  typename Container::allocator_type          allocator_type;
  typedef  typename Container::value_type              index_type;
  typedef  typename Container::reference               reference;
  typedef  typename Container::const_reference         const_reference;
  typedef  typename Container::pointer                 pointer;
  typedef  typename Container::const_pointer           const_pointer;

  typedef  typename Container::iterator                value_iterator;
  typedef  typename Container::const_iterator          value_const_iterator;
  typedef  typename Container::reverse_iterator        value_reverse_iterator;
  typedef  typename Container::const_reverse_iterator  value_const_reverse_iterator;

  typedef  std::map<K,
                    Container,
                    TypeCompare>                       value_map;

  typedef  typename value_map::key_type                key_type;
  typedef  typename value_map::value_type              value_type;
  typedef  typename value_map::iterator                type_iterator;
  typedef  typename value_map::const_iterator          type_const_iterator;
  typedef  typename value_map::reverse_iterator        type_reverse_iterator;
  typedef  typename value_map::const_reverse_iterator  type_const_reverse_iterator;

  typedef           std::pair<const K, T>              element_type;

  typedef           ValueCompare                       value_compare;
  typedef           TypeCompare                        type_compare;

public:

  class iterator : public std::iterator<std::bidirectional_iterator_tag, T, difference_type, pointer, reference>
  {
    public:
      friend class MultiplexMap<K, T, TypeCompare, ValueCompare, Container>;
      friend class self_type::const_iterator;

      iterator();
     
      bool              operator==(const iterator& x)       const;
      bool              operator!=(const iterator& x)       const;
      bool              operator==(const const_iterator& x) const;
      bool              operator!=(const const_iterator& x) const;

      iterator&         operator=(const iterator& x);
      
      iterator&         operator++();
      iterator&         operator--();
      iterator          operator++(int);
      iterator          operator--(int);
      
      typename self_type::pointer           operator->()                        const;
      typename self_type::reference         operator*()                         const;

      value_iterator    value() const { return my_value_iterator; }
      type_iterator     type()  const { return my_type_iterator;  }
    
    protected:

      iterator(MultiplexMap<K, T, TypeCompare, ValueCompare, Container>*, type_iterator, value_iterator);

      void                              advance_iterator();
      void                              retreat_iterator();
      
      type_iterator                     my_type_iterator;
      value_iterator                    my_value_iterator;

      MultiplexMap<K, T, 
                   TypeCompare, 
                   ValueCompare, 
                   Container>*          my_multiplexmap;
  };  

  class const_iterator : public std::iterator<std::bidirectional_iterator_tag, T, difference_type, pointer, const_reference>
  {
    public:
      friend class MultiplexMap<K, T, TypeCompare, ValueCompare, Container>;
      friend class self_type::iterator;

      const_iterator();
      const_iterator(const self_type::iterator& x);
     
      bool                    operator==(const self_type::iterator& x) const;
      bool                    operator!=(const self_type::iterator& x) const;
      bool                    operator==(const const_iterator& x)      const;
      bool                    operator!=(const const_iterator& x)      const;

      const_iterator&         operator=(const self_type::iterator& x);
      const_iterator&         operator=(const const_iterator& x);
      
      const_iterator&         operator++();
      const_iterator&         operator--();
      const_iterator          operator++(int);
      const_iterator          operator--(int);
      
      const_pointer           operator->()                             const;
      const_reference         operator*()                              const;

      value_const_iterator    value() const { return my_value_const_iterator; }
      type_const_iterator     type()  const { return my_type_const_iterator;  }
    
    protected:

      const_iterator(const MultiplexMap<K, T, TypeCompare, ValueCompare, Container>*, type_const_iterator, value_const_iterator);

      void                                advance_const_iterator();
      void                                retreat_const_iterator();

      type_const_iterator                 my_type_const_iterator;
      value_const_iterator                my_value_const_iterator;

      const MultiplexMap<K, T, 
                         TypeCompare, 
                         ValueCompare, 
                         Container>*      my_multiplexmap;
  };  

public:

  // Constructors.
  MultiplexMap();
  
  MultiplexMap(const type_compare& tc);

  MultiplexMap(const type_compare& tc, const value_compare& vc);

  template <class InputIterator>
  MultiplexMap(InputIterator first, InputIterator last, const type_compare& tc);

  template <class InputIterator>
  MultiplexMap(InputIterator first, InputIterator last, const type_compare& tc, const value_compare& vc);

  MultiplexMap(const MultiplexMap& mm);

  MultiplexMap&                    operator=(const MultiplexMap& mm);

  // Modifiers.
  std::pair<iterator,bool>         insert(const element_type& item, bool no_create = false);

  template <class IteratorType>
  void                             insert(IteratorType begin, IteratorType end);
  
  void                             erase_type(iterator position);
  void                             erase_type(iterator first, iterator last);
  size_type                        erase_type(const key_type& k);
  
  void                             erase_value(iterator position);
  void                             erase_value(iterator first, iterator last);
  size_type                        erase_value(const element_type& item);
  
  void                             swap(MultiplexMap& mm);
  void                             clear();
  
  // Capacity.
  bool                             type_empty()                                   const;
  size_type                        type_size()                                    const;
  size_type                        type_max_size()                                const;
  size_type                        type_count(const key_type& k)                  const;
  
  bool                             value_empty(const key_type& k)                 const;
  size_type                        value_size(const key_type& k)                  const;
  size_type                        value_max_size(const key_type& k)              const;
  size_type                        value_count(const element_type& item)          const;

  // Element access.
  Container&                       operator[](const key_type& k);
  const Container&                 operator[](const key_type& k)                  const;
  
  // Observers.
  type_compare                     type_comp()                                    const;
  value_compare                    value_comp()                                   const;

  type_iterator                    find(const key_type& k);
  type_const_iterator              find(const key_type& k)                        const;

  iterator                         find(const element_type& item);
  const_iterator                   find(const element_type& item)                 const;
  
  // Iterations.  
  type_iterator                    type_begin();
  type_const_iterator              type_begin()                                   const;

  type_iterator                    type_end();
  type_const_iterator              type_end()                                     const;

  type_reverse_iterator            type_rbegin();
  type_const_reverse_iterator      type_rbegin()                                  const;

  type_reverse_iterator            type_rend();
  type_const_reverse_iterator      type_rend()                                    const;

  value_iterator                   value_begin(const key_type& k);
  value_const_iterator             value_begin(const key_type& k)                 const;

  value_iterator                   value_end(const key_type& k);
  value_const_iterator             value_end(const key_type& k)                   const;

  value_reverse_iterator           value_rbegin(const key_type& k);
  value_const_reverse_iterator     value_rbegin(const key_type& k)                const;

  value_reverse_iterator           value_rend(const key_type& k);
  value_const_reverse_iterator     value_rend(const key_type& k)                  const;

  iterator                         begin();
  const_iterator                   begin()                                        const;

  iterator                         end();
  const_iterator                   end()                                          const;

  // Operations.
  std::pair<iterator,iterator>     equal_range(const key_type& k);
  std::pair<const_iterator,
            const_iterator>        equal_range(const key_type& k)                 const;

  iterator                         lower_bound(const key_type& k);
  const_iterator                   lower_bound(const key_type& k)                 const;
  
  iterator                         upper_bound(const key_type& k);
  const_iterator                   upper_bound(const key_type& k)                 const;

  bool                             operator==(const MultiplexMap& ms)             const;
  bool                             operator!=(const MultiplexMap& ms)             const;
  bool                             operator< (const MultiplexMap& ms)             const;
  bool                             operator<=(const MultiplexMap& ms)             const;
  bool                             operator> (const MultiplexMap& ms)             const;
  bool                             operator>=(const MultiplexMap& ms)             const;

private:

  value_compare                    my_value_compare;
  value_map                        my_multiplexer;

}; // end of class definition

//----------------------------------------------------------------------------
// Implementation of iterator
//----------------------------------------------------------------------------
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::iterator()
{
  my_multiplexmap   = NULL;
  my_type_iterator  = type_iterator();
  my_value_iterator = value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
::iterator(MultiplexMap<K, T, TypeCompare, ValueCompare, Container>* mm, type_iterator t, value_iterator v)
{
  my_multiplexmap   = mm;
  my_type_iterator  = t;
  my_value_iterator = v;
}
            
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator==(const iterator& x) const
{
  return my_type_iterator == x.my_type_iterator && my_value_iterator == x.my_value_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator!=(const iterator& x) const
{
  return !(*this == x);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator==(const const_iterator& x) const
{
  return my_type_iterator == x.my_type_iterator && my_value_iterator == x.my_value_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator!=(const const_iterator& x) const
{
  return !(*this == x);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator=(const iterator& x)
{
  my_multiplexmap   = x.my_multiplexmap;
  my_type_iterator  = x.my_type_iterator;
  my_value_iterator = x.my_value_iterator;
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator++()
{
  advance_iterator();
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator--()
{
  retreat_iterator();
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator++(int)
{
  iterator temp(*this);
  advance_iterator();
  return temp;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator--(int)
{
  iterator temp(*this);
  retreat_iterator();
  return temp;
}
      
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::pointer
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator->() const
{
  return &( operator*() );
}
      
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::reference
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::operator*() const
{
  return *my_value_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::advance_iterator()
{
  ++my_value_iterator;

  if( my_value_iterator != my_type_iterator->second.end() )
    return;

  type_iterator type_end = my_multiplexmap->type_end();

  if( my_type_iterator != type_end )
    ++my_type_iterator;

  while( my_type_iterator != type_end && !my_type_iterator->second.size() )
    ++my_type_iterator;
  
  if( my_type_iterator == type_end )
    (*this) = my_multiplexmap->end();
  else  
    my_value_iterator = my_type_iterator->second.begin();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator::retreat_iterator()
{
  if( my_value_iterator != my_type_iterator->second.begin() )
  {
    --my_value_iterator;
    return;
  }
  
  type_iterator type_begin = my_multiplexmap->type_begin();

  if( my_type_iterator != type_begin )
    --my_type_iterator;

  while( my_type_iterator != type_begin && !my_type_iterator->second.size() )
    --my_type_iterator;
  
  if( !my_type_iterator->second.size())
    (*this) = my_multiplexmap->begin();
  else  
    my_value_iterator = --my_type_iterator->second.end();
}

//----------------------------------------------------------------------------
// Implementation of const_iterator
//----------------------------------------------------------------------------

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::const_iterator()
{
  my_multiplexmap         = NULL;
  my_type_const_iterator  = type_const_iterator();
  my_value_const_iterator = value_const_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::
const_iterator::const_iterator(const MultiplexMap<K, T, TypeCompare, ValueCompare, Container>* mm, type_const_iterator t, value_const_iterator v)
{
  my_multiplexmap         = mm;
  my_type_const_iterator  = t;
  my_value_const_iterator = v;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::const_iterator(const self_type::iterator& x)
{
  my_multiplexmap         = x.my_multiplexmap;
  my_type_const_iterator  = x.my_type_iterator;
  my_value_const_iterator = x.my_value_iterator;
}
            
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator==(const self_type::iterator& x) const
{
  return my_type_const_iterator == x.my_type_const_iterator && my_value_const_iterator == x.my_value_const_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator!=(const self_type::iterator& x) const
{
  return !(*this == x);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator==(const const_iterator& x) const
{
  return my_type_const_iterator == x.my_type_const_iterator && my_value_const_iterator == x.my_value_const_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator!=(const const_iterator& x) const
{
  return !(*this == x);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator=(const self_type::iterator& x)
{
  my_multiplexmap   = x.my_multiplexmap;
  my_type_const_iterator  = x.my_type_iterator;
  my_value_const_iterator = x.my_value_iterator;
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator=(const const_iterator& x)
{
  my_multiplexmap         = x.my_multiplexmap;
  my_type_const_iterator  = x.my_type_const_iterator;
  my_value_const_iterator = x.my_value_const_iterator;
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator++()
{
  advance_const_iterator();
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator--()
{
  retreat_const_iterator();
  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator++(int)
{
  const_iterator temp(*this);
  advance_const_iterator();
  return temp;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator--(int)
{
  const_iterator temp(*this);
  retreat_const_iterator();
  return temp;
}
      
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_pointer
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator->() const
{
  return &( operator*() );
}
      
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_reference
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::operator*() const
{
  return *my_value_const_iterator;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::advance_const_iterator()
{
  ++my_value_const_iterator;

  if( my_value_const_iterator != my_type_const_iterator->second.end() )
    return;

  type_const_iterator type_end = my_multiplexmap->type_end();

  if( my_type_const_iterator != type_end )
    ++my_type_const_iterator;

  while( my_type_const_iterator != type_end && !my_type_const_iterator->second.size() )
    ++my_type_const_iterator;
  
  if( my_type_const_iterator == type_end )
    (*this) = my_multiplexmap->end();
  else  
    my_value_const_iterator = my_type_const_iterator->second.begin();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator::retreat_const_iterator()
{
  if( my_value_const_iterator != my_type_const_iterator->second.begin() )
  {
    --my_value_const_iterator;
    return;
  }
  
  type_const_iterator type_begin = my_multiplexmap->type_begin();

  if( my_type_const_iterator != type_begin )
    --my_type_const_iterator;

  while( my_type_const_iterator != type_begin && !my_type_const_iterator->second.size() )
    --my_type_const_iterator;
  
  if( !my_type_const_iterator->second.size())
    (*this) = my_multiplexmap->begin();
  else  
    my_value_const_iterator = --my_type_const_iterator->second.end();
}

//----------------------------------------------------------------------------
// Implementation of MultiplexMap
//----------------------------------------------------------------------------

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap()
{ }

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap(const type_compare& tc) : my_multiplexer(tc)
{ }

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap(const type_compare& tc, const value_compare& vc)
              : my_multiplexer(tc), my_value_compare(vc)
{ }

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
template <class InputIterator>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap(InputIterator first, InputIterator last, const type_compare& tc)
              : my_multiplexer(first, last, tc)
{ }
              
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
template <class InputIterator>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap(InputIterator first, InputIterator last, const type_compare& tc, const value_compare& vc)
              : my_multiplexer(first, last, tc), my_value_compare(vc)
{ }

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::MultiplexMap(const MultiplexMap& mm)
{
  my_multiplexer   = mm.my_multiplexer;
  my_value_compare = mm.my_value_compare;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator=(const MultiplexMap& mm)
{
  my_multiplexer   = mm.my_multiplexer;
  my_value_compare = mm.my_value_compare;

  return *this;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
std::pair<typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator,bool>
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>
::insert(const element_type& item, bool no_create)
{
  type_iterator  s = find(item.first);

  if( s == type_end() )
  {
    if(no_create)
      return std::make_pair(end(),false);
      
    Container c(my_value_compare);
    std::pair<typename Container::iterator, bool> vib = c.insert(item.second);
    typename value_map::value_type vp(item.first, c);
    std::pair<typename value_map::iterator, bool> mib = my_multiplexer.insert(vp);

    return std::make_pair( iterator(this, mib.first, vib.first), vib.second );
  }
  else
  {
    std::pair<typename Container::iterator, bool> vib = s->second.insert(item.second);

    return std::make_pair( iterator(this, s, vib.first), vib.second );
  }
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
template <class InputIterator>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::insert(InputIterator begin, InputIterator end)
{
  while( begin != end )
    insert(*begin++);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_type(iterator position)
{
  type_iterator t = position.type();
  my_multiplexer.erase(t);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_type(iterator first, iterator last)
{
  type_iterator f = first.type();
  type_iterator l = last.type();
  my_multiplexer.erase(f, l);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_type(const key_type& k)
{
  return my_multiplexer.erase(k);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_value(iterator position)
{
  type_iterator  t = position.type();
  value_iterator v = position.value();

  t->second.erase(v);  
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_value(iterator first, iterator last)
{
  type_iterator  t_first = first.type();
  value_iterator v_first = first.value();
  
  type_iterator  t_last  = last.type();
  value_iterator v_last  = last.value();
  
  if( t_first == t_last )
  {
    t_first->second.erase(v_first, v_last);
    
    return;
  }
  
  value_iterator v_end   = t_first->second.end();
  t_first->second.erase(v_first, v_end);
  
  ++t_first;
  for( type_iterator temp = t_first; temp != t_last; ++temp )
    temp->second.erase(temp->second.begin(), temp->second.end());
    
  value_iterator v_begin = t_last->second.begin();
  t_last->second.erase(v_begin, v_last);
  
  return;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::erase_value(const element_type& item)
{
  type_iterator t = find(item.first);
  if( t == type_end() )
    return 0;
  return t->second.erase(item.second);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::swap(MultiplexMap& mm)
{
  my_multiplexer.swap(mm.my_multiplexer);
  std::swap(my_value_compare, mm.my_value_compare);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
void
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::clear()
{
  my_multiplexer.clear();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_compare
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_comp() const
{
  return my_multiplexer.key_comp();
}
  
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_compare
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_comp() const
{
  return my_value_compare;
}      

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
Container&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator[](const key_type& k)
{
  return my_multiplexer[k];
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
const Container&
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator[](const key_type& k) const
{
  return my_multiplexer[k];
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_empty() const
{
  return my_multiplexer.empty();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_size() const
{
  return my_multiplexer.size();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_max_size() const
{
  return my_multiplexer.max_size();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_count(const key_type& k) const
{
  return my_multiplexer.count(k);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_empty(const key_type& k) const
{
  type_const_iterator t = find(k);

  if( t == type_end() )
    return true;
    
  return t->second.empty();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_size(const key_type& k) const
{
  type_const_iterator t = find(k);

  if( t == type_end() )
    return 0;
    
  return t->second.size();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_max_size(const key_type& k) const
{
  type_const_iterator t = find(k);

  if( t == type_end() )
    return 0;
    
  return t->second.max_size();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::size_type
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_count(const element_type& item) const
{
  type_const_iterator t = find(item.first);

  if( t == type_end() )
    return 0;
    
  return t->second.count(item.second);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
std::pair<typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator,
          typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator>
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::equal_range(const key_type& k)
{
  pair<type_iterator, type_iterator>  t = my_multiplexer.equal_range(k);

  value_iterator v_begin = t.first->second.begin();
  value_iterator v_end   = t.first->second.end();
  
  if( t.second != type_end() )
  {
    v_end = t.second->second.begin();
    return make_pair(iterator(this, t.first, v_begin), iterator(this, t.second, v_end));
  }
  
  return make_pair(iterator(this, t.first, v_begin), iterator(this, t.first, v_end));
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
std::pair<typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator,
          typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator>
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::equal_range(const key_type& k) const
{
  pair<type_const_iterator, type_const_iterator>  t = my_multiplexer.equal_range(k);

  value_const_iterator v_begin = this->value_begin(k);
  value_const_iterator v_end   = this->value_end(k);

  if( t.second != type_end() )
  {
    v_end = t.second->second.begin();
    return make_pair(const_iterator(this, t.first, v_begin), iterator(this, t.second, v_end));
  }
  
  return make_pair(const_iterator(this, t.first, v_begin), iterator(this, t.first, v_end));
}


template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::lower_bound(const key_type& k)
{
  type_iterator t = my_multiplexer.lower_bound(k);
  
  if( t == type_end() )
    return end();
    
  value_iterator vi = t->second.begin();
  return iterator(this, t, vi);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::lower_bound(const key_type& k) const
{
  type_const_iterator t = my_multiplexer.lower_bound(k);
  
  if( t == type_end() )
    return end();
    
  value_const_iterator vi = t->second.begin();
  return const_iterator(this, t, vi);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::upper_bound(const key_type& k)
{
  type_iterator t = my_multiplexer.upper_bound(k);
  
  if( t == type_end() )
    return end();

    
  value_iterator vi = t->second.begin();
  return iterator(this, t, vi);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::upper_bound(const key_type& k) const
{
  type_const_iterator t = my_multiplexer.upper_bound(k);
  
  if( t == type_end() )
    return end();
    
  value_const_iterator vi = t->second.begin();
  return const_iterator(this, t, vi);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::find(const key_type& k)
{
  return my_multiplexer.find(k);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::find(const key_type& k) const
{
  return my_multiplexer.find(k);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::find(const element_type& item)
{
  type_iterator t = find(item.first);
  
  if( t != type_end() )
  {
    value_iterator v = t->second.find(item.second);
    if( v != t->second.end() )
      return iterator(this, t, v);
  }
  
  return end();
}  

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::find(const element_type& item) const
{
  type_const_iterator t = find(item.first);
  
  if( t != type_end() )
  {
    value_const_iterator v = t->second.find(item.second);
    if( v != t->second.end() )
      return const_iterator(this, t, v);
  }
  
  return end();
}  

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_begin()
{
  return my_multiplexer.begin();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_begin()   const
{
  return my_multiplexer.begin();
}
template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_end()
{
  return my_multiplexer.end();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_end()     const
{
  return my_multiplexer.end();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_rbegin()
{
  return my_multiplexer.rbegin();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_const_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_rbegin()  const
{
  return my_multiplexer.rbegin();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_rend()
{
  return my_multiplexer.rend();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_const_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::type_rend()  const
{
  return my_multiplexer.rend();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_begin(const key_type& v)
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->begin();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_begin(const key_type& v) const
{
  type_const_iterator s = find(v);
  if( s != type_end() )
    return s->begin();
  return value_const_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_end(const key_type& v)
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->end();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_end(const key_type& v) const
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->end();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_rbegin(const key_type& v)
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->rbegin();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_const_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_rbegin(const key_type& v) const
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->rbegin();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_rend(const key_type& v)
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->rend();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_const_reverse_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::value_rend(const key_type& v) const
{
  type_iterator s = find(v);
  if( s != type_end() )
    return s->rend();
  return value_iterator();
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::begin()
{
  if( !type_size() )
    return iterator(this, type_end(), value_iterator() );     

  type_iterator t  = this->type_begin();
  type_iterator te = this->type_end();
  
  for( ; t != te && !t->second.size() ; ++t);

  if(t == te)
    return iterator(this, t, value_iterator());     
  else
    return iterator(this, t, t->second.begin() );     
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::begin() const
{
  if( !type_size() )
    return const_iterator(this, type_end(), value_const_iterator() );     
  
  type_const_iterator t  = this->type_begin();
  type_const_iterator te = this->type_end();
  
  for( ; t != te && !t->second.size() ; ++t);

  if(t == te)
    return const_iterator(this, t, value_const_iterator());     
  else
    return const_iterator(this, t, t->second.begin() );     
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::end()
{
  if( !type_size() )
    return iterator(this, type_end(), value_iterator() );     

  type_iterator tb = this->type_begin();
  type_iterator t  = --this->type_end();
  
  for( ; t != tb && !t->second.size() ; --t);

  if( !t->second.size() )
    return iterator(this, t, value_iterator());     
  else
    return iterator(this, t, t->second.end() );     
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
typename MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::const_iterator
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::end() const
{
  if( !type_size() )
    return const_iterator(this, type_end(), value_iterator() );     

  type_const_iterator tb = this->type_begin();
  type_const_iterator t  = --this->type_end();
  
  for( ; t != tb && !t->second.size() ; --t);

  if( !t->second.size() )
    return const_iterator(this, t, value_const_iterator());     
  else
    return const_iterator(this, t, t->second.end() );     
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator==(const MultiplexMap& mm) const
{
  return my_multiplexer == mm.my_multiplexer;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator!=(const MultiplexMap& mm) const
{
  return !(*this == mm);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator< (const MultiplexMap& mm) const
{
  return my_multiplexer < mm.my_multiplexer;
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator> (const MultiplexMap& mm) const
{
  return (mm < *this);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator<=(const MultiplexMap& mm) const
{
  return !(mm < *this);
}

template <class K, class T, class TypeCompare, class ValueCompare, class Container>
inline
bool
MultiplexMap<K, T, TypeCompare, ValueCompare, Container>::operator>=(const MultiplexMap& mm) const
{
  return !(*this < mm);
}

} // end of namespace SAGE

#endif
