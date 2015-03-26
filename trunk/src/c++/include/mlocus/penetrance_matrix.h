#ifndef PENETRANCE_MATRIX_H
#define PENETRANCE_MATRIX_H

//============================================================================
//  File:       sparse.h
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include <globals/config.h>
#include <assert.h>
#include <algorithm>

namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  CLASS:  penetrance_matrix
//============================================================================
//
template <class RowIndex, class ColIndex>
class penetrance_matrix
{
  protected:
  
    struct location
    {
      RowIndex row;
      ColIndex col;
      
      location() : row(0), col(0) {}
      location(RowIndex r, ColIndex c) : row(r), col(c) {}
      location(const location& l) : row(l.row), col(l.col) {}

      location& operator=(const location& l)
      {
        row = l.row;
        col = l.col;
      }
      
      bool operator<(const location& l) const
      {
        return (row < l.row) || (row == l.row && col < l.col);
      }
    };

    typedef std::map<location, double> element_map;

  public:
  
    typedef typename element_map::iterator       row_iterator;

  protected:
  
    typedef typename element_map::value_type        map_value;
    typedef std::pair<row_iterator, row_iterator>   iterator_pair;
    typedef std::vector<iterator_pair>              col_vector;

  public:
    penetrance_matrix();
    penetrance_matrix(RowIndex row_size);
    penetrance_matrix(const penetrance_matrix<RowIndex, ColIndex>&);
    virtual ~penetrance_matrix();

    penetrance_matrix&  operator =(const penetrance_matrix<RowIndex, ColIndex>&);

    row_iterator    row_begin(RowIndex row) const;
    row_iterator    row_end(RowIndex row) const;

    uint    row_elements(RowIndex r) const;

    double       operator ()(RowIndex row, ColIndex col) const;
    double       default_value() const;
   
    penetrance_matrix&  clear();
    penetrance_matrix&  clear_row(RowIndex row);
    penetrance_matrix&  resize(RowIndex row_size);
    penetrance_matrix&  set(RowIndex row, ColIndex col, const double& val);
    penetrance_matrix&  remove(RowIndex row, ColIndex col);
    penetrance_matrix&  set_default_value(const double& val);
    penetrance_matrix&  swap(penetrance_matrix<RowIndex, ColIndex>& M);

    // Returns if the value is set to anything.

    bool set(RowIndex row, ColIndex col) const;

  private:
    double      my_default;
    col_vector  my_data;

    element_map  my_elements;
    element_map  my_set_default;
};

//============================================================================
//  IMPLEMENTATION: penetrance_matrix
//============================================================================
//
template <class RowIndex, class ColIndex>
penetrance_matrix<RowIndex, ColIndex>::penetrance_matrix()
  : my_default(), my_data(1)
{
  my_data[0] = make_pair(my_elements.end(), my_set_default.end());
}

template <class RowIndex, class ColIndex>
penetrance_matrix<RowIndex, ColIndex>::penetrance_matrix(RowIndex row_size)
  : my_default(), my_data()
{
    my_data.resize(row_size+1, make_pair(my_elements.end(), my_set_default.end()));
}

template <class RowIndex, class ColIndex>
penetrance_matrix<RowIndex, ColIndex>::penetrance_matrix(const penetrance_matrix<RowIndex, ColIndex>& M)
  : my_default(M.my_default), 
    my_data(M.my_data), my_elements(M.my_elements)
{
    uint r = 0;

    for( ; r < my_data.size()-1; ++r)
    {
        // This was fixed 2002-05-21.  See SCR #262 in Mantis for details

        if(M.my_data[r].first == M.my_elements.end())
        {
          row_iterator i = my_elements.end();
          row_iterator j = my_set_default.end();
          
          my_data[r] = make_pair(i,j);
        }
        else
        {
          row_iterator i = my_elements.find(my_data[r].first->first);
          row_iterator j = my_set_default.find(my_data[r].second->first);
          
          my_data[r] = make_pair(i,j);
        }
    }

    my_data.back() = make_pair(my_elements.end(), my_set_default.end());
}

template <class RowIndex, class ColIndex>
penetrance_matrix<RowIndex, ColIndex>::~penetrance_matrix() 
{
}

//----------
//
template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::operator =(const penetrance_matrix<RowIndex, ColIndex>& M)
{
    if (&M != this)
    {
        penetrance_matrix<RowIndex, ColIndex> tmp(M);

        swap(tmp);
    }
    return *this;
}

//----------
//
template <class RowIndex, class ColIndex> double
penetrance_matrix<RowIndex, ColIndex>::operator ()(RowIndex row, ColIndex col) const
{
    typename element_map::const_iterator rf = my_elements.find(location(row, col));

    if (rf != my_elements.end())
    {
        return rf->second;
    }

    return my_default;
}

template <class RowIndex, class ColIndex> inline double
penetrance_matrix<RowIndex, ColIndex>::default_value() const
{
    return my_default;
}

template <class RowIndex, class ColIndex> inline uint
penetrance_matrix<RowIndex, ColIndex>::row_elements(RowIndex row) const
{
  uint count  = 0;
  
  row_iterator i = my_data[row].first;
  
  for( ; i != my_data[row+1].first; ++count, ++i);
  
  return count;
}

//----------
//
template <class RowIndex, class ColIndex> inline typename penetrance_matrix<RowIndex, ColIndex>::row_iterator
penetrance_matrix<RowIndex, ColIndex>::row_begin(RowIndex row) const
{
    return my_data[row].first;
}

template <class RowIndex, class ColIndex> inline typename penetrance_matrix<RowIndex, ColIndex>::row_iterator
penetrance_matrix<RowIndex, ColIndex>::row_end(RowIndex row) const
{
    return my_data[row+1].first;
}

//----------
//
template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::clear()
{
    for (size_t i = 0;  i < my_data.size();  ++i)
    {
        my_data[i] = make_pair(my_elements.end(), my_set_default.end());
    }

    my_elements.clear();
    my_set_default.clear();

    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::clear_row(RowIndex row)
{
    my_elements.erase(my_data[row].first, my_data[row+1].first);
    my_set_default.erase(my_data[row].second, my_data[row+1].second);

    RowIndex r = row;

    row_iterator i = my_data[row].first;
    
    while(r >= 0 && my_data[r].first == i)
    {
      my_data[r].first = my_data[row+1].first;
      --r;
    }

    r = row;

    while(r >= 0 && my_data[r].second == i)
    {
        my_data[r].second = my_data[row+1].second;
        --r;
    }
    
    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::resize(RowIndex row_size)
{
    my_data.resize(row_size+1, make_pair(my_elements.end(), my_set_default.end()));

    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::set(RowIndex row, ColIndex col, const double& val)
{
    if(row < 0 || row >= (int) my_data.size()-1)
      return *this;

    location l = location(row, col);

    if(val != my_default)
    {
      my_elements[l] = val;

      row_iterator i = my_elements.find(l);

      while(row >= 0 && (my_data[row].first == my_elements.end() || l < my_data[row].first->first))
      {
        my_data[row].first = i;
        --row;
      }
      row_iterator j = my_set_default.find(l);

      if(j != my_set_default.end()) my_set_default.erase(j);
    }
    else
    {
      my_set_default[l] = val;

      row_iterator i = my_set_default.find(l);

      while(row >= 0 && (my_data[row].second == my_set_default.end() || l < my_data[row].second->first))
      {
        my_data[row].second = i;
        --row;
      }
      row_iterator j = my_elements.find(l);

      if(j != my_elements.end()) my_elements.erase(j);
    }      

    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::remove(RowIndex row, ColIndex col)
{
    location l = location(row, col);

    row_iterator i = my_elements.find(l);
    
    if(i != my_elements.end())
    {
      while(row >= 0 && my_data[row].first == i)
      {
        ++my_data[row].first;
        --row;
      }

      my_elements.erase(i);

      return *this;
    }
    
    row_iterator j = my_set_default.find(l);
    
    if(j != my_set_default.end())
    {
      while(row >= 0 && my_data[row].second == j)
      {
        ++my_data[row].second;
        --row;
      }

      my_set_default.erase(j);
    }
    
    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::set_default_value(const double& val)
{
    element_map tmp;

    // Find all the elements that are the new default value and create the set
    {
      row_iterator b = my_elements.begin();
      row_iterator e = my_elements.end();
    
      for( ; b != e; )
      {
        if(b->second == val)
        {
          location l = b->first;

          tmp[l] = b->second;

          remove(l.row, l.col);
          
          b = my_elements.lower_bound(l);
        }
        else
          ++b;
      }
    }

    // Insert the elements that were the old default values into the map
    {
      row_iterator b = my_set_default.begin();
      row_iterator e = my_set_default.end();
      
      for( ; b != e; ++b)
      {
        my_elements[b->first] = b->second;
      }
    }

    my_set_default.swap(tmp);

    my_default = val;

    return *this;
}

template <class RowIndex, class ColIndex> penetrance_matrix<RowIndex, ColIndex>&
penetrance_matrix<RowIndex, ColIndex>::swap(penetrance_matrix<RowIndex, ColIndex>& M)
{
    std::swap(my_default, M.my_default);
    my_data.swap(M.my_data);
    my_elements.swap(M.my_elements);
    my_set_default.swap(M.my_set_default);

    return *this;
}

template <class RowIndex, class ColIndex> bool
penetrance_matrix<RowIndex, ColIndex>::set(RowIndex row, ColIndex col) const
{
  location l(row, col);

  return my_elements.find(l)    != my_elements.end() || 
         my_set_default.find(l) != my_set_default.end();
}

} // End namespace MLOCUS
} // End namespace SAGE

#endif
