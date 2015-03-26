#ifndef MARKED_VECTOR
#define MARKED_VECTOR

#include <memory>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <iostream>
using namespace std;

namespace SAGE
{

template <class T, class Allocator = std::allocator<T> >
class marked_vector
{
public:
    typedef std::vector<T,Allocator>                storage_vector;
    
    typedef storage_vector::pointer                 pointer;
    typedef storage_vector::const_pointer           const_pointer;
    typedef storage_vector::reference               reference;
    typedef storage_vector::const_reference         const_reference;    
    typedef storage_vector::size_type               size_type;
    typedef storage_vector::difference_type         difference_type;
    typedef storage_vector::value_type              value_type;
    typedef storage_vector::allocator_type          allocator_type;
        
    typedef storage_vector::iterator                iterator;
    typedef storage_vector::const_iterator          const_iterator;
    typedef storage_vector::reverse_iterator        reverse_iterator;
    typedef storage_vector::const_reverse_iterator  const_reverse_iterator;

    typedef std::vector<size_t>                     mark_vector;
    typedef mark_vector::iterator                   mark_iterator;
    typedef mark_vector::const_iterator             mark_const_iterator;
    typedef mark_vector::reverse_iterator           mark_reverse_iterator;
    typedef mark_vector::const_reverse_iterator     mark_const_reverse_iterator;

     explicit marked_vector(const Allocator& alloc = Allocator()) 
          : my_storage(alloc) { }

     explicit marked_vector(size_type n)
          : my_storage(n) { }

     explicit marked_vector(size_type n, const T& value, 
                              const Allocator& alloc = Allocator())
          : my_storage(n, value, alloc) { }

     marked_vector(const std::vector<T, Allocator>& x) 
          : my_storage(x) { }

     template<class InputIterator>
     marked_vector(InputIterator first, InputIterator last, 
                         const Allocator& alloc = Allocator())
          : my_storage(first, last, alloc) { }

     template<class InputIterator>
     void assign(InputIterator first, InputIterator last)
     {
       my_marks.clear();
       my_storage.assign(first,last);
     }

     // This is what the standard draft says
     template <class Size, class Type>
     void assign(Size n, const Type& t = Type())
     {
       my_marks.clear();
       my_storage.assign(n, t);
     }

     allocator_type get_allocator() const 
     { 
       return my_storage.get_allocator(); 
     }

     ~marked_vector() { }

     marked_vector<T, Allocator>& operator=(const std::vector<T, Allocator>& x)
     {
       my_marks.clear();
       my_storage = x;
     }

     marked_vector<T, Allocator>& operator=(const marked_vector<T, Allocator>& x)
     {
       my_marks   = x.my_marks;
       my_storage = x.my_storage;
     }
     
     iterator                     begin()              { return my_storage.begin();  }
     const_iterator               begin()       const  { return my_storage.begin();  }
     iterator                     end()                { return my_storage.end();    }
     const_iterator               end()         const  { return my_storage.end();    }
     reverse_iterator             rbegin()             { return my_storage.rbegin(); }
     const_reverse_iterator       rbegin()      const  { return my_storage.rbegin(); }
     reverse_iterator             rend()               { return my_storage.rend();   }
     const_reverse_iterator       rend()        const  { return my_storage.rend();   } 

     mark_iterator                mark_begin()         { return my_marks.begin();  }
     mark_const_iterator          mark_begin()  const  { return my_marks.begin();  }
     mark_iterator                mark_end()           { return my_marks.end();    }
     mark_const_iterator          mark_end()    const  { return my_marks.end();    }
     mark_reverse_iterator        mark_rbegin()        { return my_marks.rbegin(); }
     mark_const_reverse_iterator  mark_rbegin() const  { return my_marks.rbegin(); }
     mark_reverse_iterator        mark_rend()          { return my_marks.rend();   }
     mark_const_reverse_iterator  mark_rend()   const  { return my_marks.rend();   }

     size_type size()       const  { return my_storage.size(); }
     size_type mark_size()  const  { return my_marks.size();   }
     size_type max_size()   const  { return std::min(my_storage.max_size(),
                                                     my_marks.max_size()); }

     size_type mark_size(const_iterator position) const
     {
       d = std::distance(my_storage.begin(), position);

       mark_iterator m = upper_bound( my_marks.begin(),
                                      my_marks.end(),
                                      d );
     
       size_type e = size();
       
       if( m != my_marks.end() )
         e = *m;
         
       return e-d;
     }

     size_type mark_size(mark_const_iterator position) const
     {
       if(position == mark_end())
         return 0;

       if( (position+1) == mark_end())
         return size() - *position;

       return *(position+1) - *position;
     }
       
     void resize(size_type new_size) 
     { 
       resize(new_size, T());
     }

     void resize(size_type new_size, T c)
     {
       if(new_size > my_storage.size())
         my_storage.resize(new_size, c);
       else
         erase( my_storage.begin()+new_size, my_storage.end() );
     }

     size_type capacity()      const { return my_storage.capacity(); }
     size_type mark_capacity() const { return my_marks.capacity();   }
     
     bool empty()              const { return my_storage.empty();    }

     void reserve(size_type n)       { my_storage.reserve(n);        }
     void mark_reserve(size_type n)  { my_marks.reserve(n);          }

    // These should not be range checked.
     reference       operator[] (size_type n)       { return my_storage[n]; }
     const_reference operator[] (size_type n) const { return my_storage[n]; }

     const_reference at(size_type n)          const { return my_storage.at(n); }
     reference       at(size_type n)                { return my_storage.at(n); }

     reference         front()        { return my_storage.front(); }
     const_reference   front() const  { return my_storage.front(); }
     reference         back()         { return my_storage.back();  }
     const_reference   back()  const  { return my_storage.back();  }

     void push_back(const T& x)      { my_storage.push_back(x);   }

     void mark_back()                
     { 
       if( !my_storage.size() )
         return;

       // Check if we are already marked
       if( my_marks.size() && my_marks.back() == my_storage.size() - 1 )
         return;

       my_marks.push_back( my_storage.size() - 1 );
     }

     void pop_back()
     {
       if( my_marks.size() && my_storage.size() &&
           my_marks.back() == (my_storage.end() - 1) )
           my_marks.pop_back();

       my_storage.pop_back();
     }

     
     void mark(iterator position)
     {
       mark(std::distance(my_storage.begin(), position));
     }

     void mark(size_type d)
     {
       mark_iterator m = lower_bound( my_marks.begin(),
                                      my_marks.end(),
                                      d );
       if( m == my_marks.end() )
         my_marks.push_back(d);
       else if( *m != d )
         my_marks.insert(m, d);
     //else do nothing if *m == d
     }

     void unmark(iterator position)
     {
       unmark(std::distance(my_storage.begin(), position));
     }

     void unmark(size_type d)
     {
       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );

       if( m != my_marks.end() && *m == d )
         my_marks.erase(m);
     }

     void flip_mark(iterator position)
     {
       flip_mark(std::distance(my_storage.begin(), position));
     }

     void flip_mark(size_type d)
     {
       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );
       if( m == my_marks.end() )     
         my_marks.push_back(d);
       else if( *m != d )
         my_marks.insert(m, d);
       else
         my_marks.erase(m);
     }

     iterator insert(iterator position) 
     { 
       return insert(position,T()); 
     }

     iterator insert(iterator position, const T& x)
     {
       if( position == my_storage.end() )
       {
         push_back(x);
         return my_storage.end() - 1;
       }  

       const size_type d = std::distance(my_storage.begin(), position);

       iterator rv = my_storage.insert(position, x);

       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );

       if( m == my_marks.end() )
         return rv;
                                         
       const size_type b = std::distance(my_marks.begin(), m);
       const size_type e = my_marks.size();
       for( size_type i = b; i < e; ++i )
         ++my_marks[i];

       return rv;
     }

     template <class InputIterator>
     void insert(iterator position, 
                  InputIterator first, InputIterator last)
     {
       if( position == my_storage.end() )
       {
         my_storage.insert(position, first, last);
         return;
       }  

       const difference_type d = distance(my_storage.begin(), position);
       const size_type old_size = my_storage.size();
             
       my_storage.insert(position, first, last);

       const difference_type new_elements = my_storage.size() - old_size; 

       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );

       if( m == my_marks.end() )
         return;
                                         
       const size_type b = std::distance(my_marks.begin(), m);
       const size_type e = my_marks.size();
       for( size_type i = b; i < e; ++i )
         my_marks[i] += new_elements;
     }

     void insert(iterator position, size_type n, const T& x)
     {
       if( position == my_storage.end() )
       {
         my_storage.insert(position, first, last);
         return;
       }  

       const difference_type d = distance(my_storage.begin(), position);
             
       my_storage.insert(position, n, x);

       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );

       if( m == my_marks.end() )
         return;
                                         
       const size_type b = std::distance(my_marks.begin(), m);
       const size_type e = my_marks.size();
       for( size_type i = b; i < e; ++i )
         my_marks[i] += n;
     }

     iterator erase(iterator position, bool remove_mark = false)
     {
       const size_type d = distance(my_storage.begin(), position);

       iterator rv = my_storage.erase(position);

       mark_iterator m = lower_bound( my_marks.begin(), 
                                      my_marks.end(), 
                                      d );

       if( m == my_marks.end() )
         return rv;

       if(*m == d) // We are erasing a position with a mark
       {
         // Check to see if we are required to remove a mark, or if we can't
         // move the mark since the next mark would be a duplicate.
         if(remove_mark || *m == my_storage.size() - 1 ||
           (m+1 != my_marks.end() && *(m+1) == d+1 ))
           m = my_marks.erase(m);
         else         // else we do do not change *m below
           ++m;
       }
                                         
       const size_type b = std::distance(my_marks.begin(), m);
       const size_type e = my_marks.size();
       for( size_type i = b; i < e; ++i )
         --my_marks[i];

       return rv;
     }

     iterator erase(iterator first, iterator last, bool remove_mark = false)
     {
       // Handle the case where nothing is erased
       if(first == last)
         return first;

       const size_type df = distance(my_storage.begin(), first);
       const size_type dl = distance(my_storage.begin(), last);

       const size_type old_size = my_storage.size();

       iterator rv = my_storage.erase(first, last);

       const difference_type deleted_elements = old_size - my_storage.size();

       mark_iterator mb = lower_bound( my_marks.begin(), 
                                       my_marks.end(), 
                                       df );

       mark_iterator me = lower_bound( mb, 
                                       my_marks.end(), 
                                       dl-1 );

       // Return if there are no marks to adjust
       if( mb == my_marks.end() )
         return rv;
       
       // We are erasing a range with marks in it
       if(mb != me) 
       {
         // Check to see if we are required to remove a mark or if
         // we can't bump the mark since the next mark would be a duplicate 
         if( me == my_marks.end() )
           me = my_marks.erase(mb,me);
         else if(remove_mark || dl >= old_size || *(me+1) == dl )
           me = my_marks.erase(mb,me+1);
         else
         {
           // erase all but the last marker & move it to the end of the region
           me = my_marks.erase(mb, me);
           *me = dl;
         }
       }

       const size_type b = std::distance(my_marks.begin(), me);
       const size_type e = my_marks.size();
       for( size_type i = b; i < e; ++i )
         my_marks[i] -= deleted_elements;

       return rv;
     }

     void swap(std::vector<T, Allocator>& x)
     {
          my_marks.clear();
          my_storage.swap(x);
     }

     void swap(marked_vector<T, Allocator>& x)
     {
          my_marks.swap(x.my_marks);
          my_storage.swap(x.my_storage);
     }

     void clear()  
     {
         my_marks.clear();
         my_storage.clear();
     }

private:
    storage_vector     my_storage;
    mark_vector        my_marks;
};

template <class T, class Allocator>
inline bool 
operator== (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{
     return x.size() == y.size() &&
            equal (x.begin(), x.end(), y.begin());
}

template <class T, class Allocator>
inline bool 
operator< (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{
     return lexicographical_compare (x.begin(), x.end(), 
                                     y.begin(), y.end());
}

template <class T, class Allocator>
inline bool 
operator!= (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (x == y); }
 
template <class T, class Allocator>
inline bool 
operator>= (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (x < y); }
 
template <class T, class Allocator>
inline bool 
operator> (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return (y < x); }
 
template <class T, class Allocator>
inline bool 
operator<= (const std::vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (y < x); }

template <class T, class Allocator>
inline bool 
operator== (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{
     return x.size() == y.size() &&
            equal (x.begin(), x.end(), y.begin());
}

template <class T, class Allocator>
inline bool 
operator< (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{
     return lexicographical_compare (x.begin(), x.end(), 
                                     y.begin(), y.end());
}

template <class T, class Allocator>
inline bool 
operator!= (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{ return ! (x == y); }
 
template <class T, class Allocator>
inline bool 
operator>= (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{ return ! (x < y); }
 
template <class T, class Allocator>
inline bool 
operator> (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{ return (y < x); }
 
template <class T, class Allocator>
inline bool 
operator<= (const marked_vector<T, Allocator>& x, const std::vector<T, Allocator>& y)
{ return ! (y < x); }
 
template <class T, class Allocator>
inline bool 
operator== (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{
     return x.size() == y.size() &&
            equal (x.begin(), x.end(), y.begin());
}

template <class T, class Allocator>
inline bool 
operator< (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{
     return lexicographical_compare (x.begin(), x.end(), 
                                     y.begin(), y.end());
}

template <class T, class Allocator>
inline bool 
operator!= (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (x == y); }
 
template <class T, class Allocator>
inline bool 
operator>= (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (x < y); }
 
template <class T, class Allocator>
inline bool 
operator> (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return (y < x); }
 
template <class T, class Allocator>
inline bool 
operator<= (const marked_vector<T, Allocator>& x, const marked_vector<T, Allocator>& y)
{ return ! (y < x); }

} 

#endif
