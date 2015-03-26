#ifndef INHERITANCE_H
#define INHERITANCE_H

//============================================================================
//  File:       inheritance_model.h
//
//  Author:     Geoff Wedig
//
//  History:    Version 0.01
//
//  Notes:      Inheritance Models (marker models) for SAGE
//
//  Copyright (c) 2000 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "LSF/parse_ops.h"
#include "containers/indexed_map.h"
#include "mlocus/penmodel.h"

namespace SAGE   {
namespace MLOCUS {

/// The penetrance model provides all the information we currently need. 
/// Thus, we're using it as a surrogate.  Eventually, a class may be made to
/// replace it, but not yet.
typedef penetrance_model                       inheritance_model;

/// The inheritance_model_map is a doubly indexed container of inheritance
/// models. It is indexed on both the name of the inheritance model and on
/// the index (ie, numerical sequence) with which it was entered.  Insertion and
/// erasure operations are provided, though note that they invalidate iterators.  In
/// addition, erasure may reindex elements.  Be careful when using these features.

class inheritance_model_map
{
    typedef indexed_map<string, inheritance_model> inheritance_map;
    
public:
    typedef inheritance_map::size_type          size_type;
    typedef inheritance_map::key_type           key_type;
    typedef inheritance_map::value_type         value_type;
    typedef size_type                           index_type;

    typedef inheritance_map::indexed_iterator               index_iterator;
    typedef inheritance_map::indexed_const_iterator         index_const_iterator;
//    typedef inheritance_map::indexed_const_reverse_iterator index_const_reverse_iterator;
//    typedef inheritance_map::indexed_reverse_iterator       index_reverse_iterator;

    typedef inheritance_map::ordered_iterator               name_iterator;
    typedef inheritance_map::ordered_const_iterator         name_const_iterator;
//    typedef inheritance_map::ordered_const_reverse_iterator name_const_reverse_iterator;
//    typedef inheritance_map::ordered_reverse_iterator       name_reverse_iterator;

    explicit inheritance_model_map() { }

    template <class InputIterator>
    inheritance_model_map(InputIterator first, InputIterator last)
    {
	insert( first, last );
    }

    inheritance_model_map(const inheritance_model_map& x) : my_map(x.my_map)
    {}

    inheritance_model_map& operator=(const inheritance_model_map& x)
    {
      if(this == &x) return *this;

      my_map = x.my_map;
      return *this;
    }

    // name iterators:
    name_iterator                name_begin()         { return my_map.ordered_begin();  }
    name_const_iterator          name_begin()   const { return my_map.ordered_begin();  }
    name_iterator                name_end()           { return my_map.ordered_end();    }
    name_const_iterator          name_end()     const { return my_map.ordered_end();    }
//    name_reverse_iterator        name_rbegin()        { return my_map.ordered_rbegin(); }
//    name_const_reverse_iterator  name_rbegin()  const { return my_map.ordered_rbegin(); }
//    name_reverse_iterator        name_rend()          { return my_map.ordered_rend();   }
//    name_const_reverse_iterator  name_rend()    const { return my_map.ordered_rend();   }

    // index iterators:
    index_iterator                index_begin()         { return my_map.indexed_begin();  }
    index_const_iterator          index_begin()   const { return my_map.indexed_begin();  }
    index_iterator                index_end()           { return my_map.indexed_end();    }
    index_const_iterator          index_end()     const { return my_map.indexed_end();    }
//    index_reverse_iterator        index_rbegin()        { return my_map.indexed_rbegin(); }
//    index_const_reverse_iterator  index_rbegin()  const { return my_map.indexed_rbegin(); }
//    index_reverse_iterator        index_rend()          { return my_map.indexed_rend();   }
//    index_const_reverse_iterator  index_rend()    const { return my_map.indexed_rend();   }
 
    // capacity:
    bool          empty()     const { return my_map.empty(); }
    size_type     size()      const { return my_map.size();  }
    size_type     max_size()  const { return my_map.max_size(); }

    // accessors:
    const string& name(size_type s) const
    {
      return my_map.key(s);
    }
    size_type index(const string& s) const
    {
      return my_map.index(toUpper(s));
    }

    // modifiers:
    inheritance_model& operator[](const key_type& x)
    {
      inheritance_model& i = my_map[toUpper(x)];
      if(i.name().size() == 0) /*lint -e{534} */ i.set_name(x);

      return i;
    }

    inheritance_model& operator[](size_type s)
    {
      return my_map[s];
    }
    const inheritance_model& operator[](size_type s) const
    {
      return my_map[s];
    }

    void push_back (const value_type& x)
    {
      my_map.push_back(x);
    }

    void pop_back ()
    {
      my_map.pop_back();
    }
    
    std::pair<name_iterator, bool> insert(const value_type& x)
    {
      return my_map.insert(x);
    }

    name_iterator insert(const name_iterator& position, const value_type& x)
    {
      return my_map.insert(position, x);
    }
      
    template <class InputIterator> 
    void insert(InputIterator first, InputIterator last)
    {
      my_map.insert(first,last);
    }

    void erase(const key_type &x)
    {
      my_map.erase(toUpper(x));
    }
    
    void erase(const name_iterator& position)
    {
      my_map.erase(position);
    }

    void erase(const index_iterator& position)
    {
      my_map.erase(position);
    }

/*
    void erase(name_iterator first, name_iterator last)
    {
      my_map.erase(first, last);
    }

    void erase(index_iterator first, index_iterator last)
    {
      my_map.erase(first, last);
    }
*/
    void clear()
    {
      my_map.clear();
    }
    
    // [lib.map.ops] set operations:
    index_iterator index_find(const key_type& x)
    {
      return my_map.indexed_find(toUpper(x));
    }

    index_const_iterator index_find(const key_type& x) const
    {
      return my_map.indexed_find(toUpper(x));
    }

    name_iterator name_find(const key_type& x)
    {
      return my_map.ordered_find(toUpper(x));
    }
    
    name_const_iterator name_find(const key_type& x) const
    {
      return my_map.ordered_find(toUpper(x));
    }

private:

    inheritance_map my_map;
};

// A couple of useful functions
void inheritance_model_test(ostream& o, const inheritance_model& m, const char* str);
void inheritance_model_print(ostream& o, const inheritance_model& m, const char* str);


} // End namespace MLOCUS
} // End namespace SAGE

#endif
