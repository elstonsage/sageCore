#ifndef TRIANGLEMATRIX_H
#define TRIANGLEMATRIX_H

//****************************************************************************
//* File:      trimatrix.h                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//* Reviewer:  Kevin Jacobs                                                  *
//* Version:   1.0                                                           *
//*                                                                          *
//* History:   v.1.0 - 1999/11/13 - reviewed by kbj                          *
//*                                                                          *
//* Purpose:   This header file defines a class template to represent a      *
//*            symmetric matrix.  It includes a complete matrix view as      *
//*            upper and lower traingular views.                             *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include <memory>
#include <stdexcept>

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     TriangleMatrix                                               ~
// ~                                                                         ~
// ~ Purpose:   Define triangle matrix container class.                      ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class T, class Allocator = std::allocator<T> >
class TriangleMatrix 
{
  public:
     class iterator;
     class const_iterator;
     friend class iterator;
     friend class const_iterator;

     enum matrix_shape { FULL=0, LOWER_TRIANGLE, UPPER_TRIANGLE };

     typedef std::vector<T, Allocator>                     vector_type;
     typedef typename vector_type::value_type              value_type;

     typedef typename vector_type::reference               reference;
     typedef typename vector_type::const_reference         const_reference;
     typedef typename vector_type::size_type               size_type;
     typedef typename vector_type::difference_type         difference_type;

     typedef typename vector_type::allocator_type          allocator_type;
     typedef typename vector_type::pointer                 pointer;
     typedef typename vector_type::const_pointer           const_pointer;

  protected:
     typedef typename vector_type::iterator                value_iterator;
     typedef typename vector_type::const_iterator          value_const_iterator;
     typedef typename vector_type::reverse_iterator        value_reverse_iterator;
     typedef typename vector_type::const_reverse_iterator  value_const_reverse_iterator;

     vector_type                                  my_vector;
     size_type                                    my_rank;

  public:
     class iterator_base
     {
          protected:

	       iterator_base(size_type r, size_type c);
               iterator_base(size_type r, size_type c, int r_d, int c_d, bool r_m);
	       iterator_base(size_type r, size_type c, int r_d, int c_d, bool r_m, matrix_shape shp);

	       size_type         my_row;
	       size_type         my_column;
	       int               my_row_direction;
	       int               my_column_direction;
	       bool              my_row_major;
	       matrix_shape      my_shape;

	       void              advance_iterator(size_type, bool);

	  public:

	       iterator_base();
	      ~iterator_base();

               size_type         row()                               const;
               size_type         column()                            const;
               int               row_direction()                     const;
               int               column_direction()                  const;
               bool              row_major()                         const;
	       matrix_shape      shape()                             const;

               void              set_row(size_type);
               void              set_column(size_type);
               void              set_row_direction(int);
               void              set_column_direction(int);
               void              set_row_major(bool);
	       void              set_shape(matrix_shape);
     };

     class iterator : public iterator_base
     {
          public:

               typedef iterator_base base_type;

	       friend class TriangleMatrix<T, Allocator>;	  
	       
	       iterator();
	      ~iterator();

	       bool              operator==(const iterator_base& x)       const;
               bool              operator!=(const iterator_base& x)       const;

	       iterator&         operator++();
	       iterator&         operator--();
	       iterator          operator++(int);
	       iterator          operator--(int);

	       pointer           operator->()                             const;
               reference         operator*()                              const;

          protected:

	       iterator(TriangleMatrix<T,Allocator>*, size_type r, size_type c);
               iterator(TriangleMatrix<T,Allocator>*, size_type r, size_type c, int r_d, int c_d, bool r_m);
	       iterator(TriangleMatrix<T,Allocator>*, size_type r, size_type c, int r_d, int c_d, bool r_m, matrix_shape shp);

               TriangleMatrix<T,Allocator>*                               my_matrix;
     };

     class const_iterator : public iterator_base
     {
	  public:

	       typedef iterator_base base_type;

	       friend class TriangleMatrix<T, Allocator>;	  

	       const_iterator();
	      ~const_iterator();

	       bool              operator==(const iterator_base& x)       const;
	       bool              operator!=(const iterator_base& x)       const;

	       const_iterator&   operator++();
	       const_iterator&   operator--();

	       const_iterator    operator++(int);
	       const_iterator    operator--(int);

	       const_pointer     operator->()                             const;
               const_reference   operator*()                              const;

          protected:

	       const_iterator(const TriangleMatrix<T,Allocator>*, size_type r, size_type c);
               const_iterator(const TriangleMatrix<T,Allocator>*, size_type r, size_type c, int r_d, int c_d, bool r_m);	       
               const_iterator(const TriangleMatrix<T,Allocator>*, size_type r, size_type c, int r_d, int c_d, bool r_m, matrix_shape shp);

               const TriangleMatrix<T,Allocator>*                         my_matrix;
     };

     typedef      iterator            reverse_iterator;
     typedef      const_iterator      const_reverse_iterator;

  public:

     explicit TriangleMatrix(const Allocator& alloc = Allocator());
     explicit TriangleMatrix(size_type n);
     explicit TriangleMatrix(size_type n, const T& value, const Allocator& alloc = Allocator());

     TriangleMatrix(const TriangleMatrix<T,Allocator>& x);
    ~TriangleMatrix();

     TriangleMatrix<T, Allocator>&  operator=(const TriangleMatrix<T,Allocator>& x);

     void                           assign(size_type n, const T& t);

     template<class InputIterator>
     void                           assign(InputIterator first, InputIterator last);

     iterator                       row_begin();
     const_iterator                 row_begin()                            const;

     iterator                       row_begin(size_type i);
     const_iterator                 row_begin(size_type i)                 const;

     iterator                       row_secondary_begin();
     const_iterator                 row_secondary_begin()                  const;

     iterator                       row_secondary_begin(size_type i);
     const_iterator                 row_secondary_begin(size_type i)       const;

     reverse_iterator               row_secondary_rbegin();
     const_reverse_iterator         row_secondary_rbegin()                 const;

     reverse_iterator               row_secondary_rbegin(size_type i);
     const_reverse_iterator         row_secondary_rbegin(size_type i)      const;

     reverse_iterator               row_rbegin();
     const_reverse_iterator         row_rbegin()                           const;

     reverse_iterator               row_rbegin(size_type i);
     const_reverse_iterator         row_rbegin(size_type i)                const;

     iterator                       column_begin();
     const_iterator                 column_begin()                         const;

     iterator                       column_begin(size_type j);
     const_iterator                 column_begin(size_type j)              const;

     iterator                       column_secondary_begin();
     const_iterator                 column_secondary_begin()               const;

     iterator                       column_secondary_begin(size_type j);       
     const_iterator                 column_secondary_begin(size_type j)    const;

     reverse_iterator               column_secondary_rbegin();
     const_reverse_iterator         column_secondary_rbegin()              const;

     reverse_iterator               column_secondary_rbegin(size_type j);
     const_reverse_iterator         column_secondary_rbegin(size_type j)   const;

     reverse_iterator               column_rbegin();
     const_reverse_iterator         column_rbegin()                        const;

     reverse_iterator               column_rbegin(size_type j);
     const_reverse_iterator         column_rbegin(size_type j)             const;

     iterator                       row_end();
     const_iterator                 row_end()                              const;

     iterator                       row_end(size_type i);
     const_iterator                 row_end(size_type i)                   const;

     iterator                       row_secondary_end();
     const_iterator                 row_secondary_end()                    const;

     iterator                       row_secondary_end(size_type i);
     const_iterator                 row_secondary_end(size_type i)         const;

     reverse_iterator               row_secondary_rend();
     const_reverse_iterator         row_secondary_rend()                   const;

     reverse_iterator               row_secondary_rend(size_type i);
     const_reverse_iterator         row_secondary_rend(size_type i)        const;

     reverse_iterator               row_rend();
     const_reverse_iterator         row_rend()                             const;

     reverse_iterator               row_rend(size_type i);
     const_reverse_iterator         row_rend(size_type i)                  const;

     iterator                       column_end();
     const_iterator                 column_end()                           const;

     iterator                       column_end(size_type j);
     const_iterator                 column_end(size_type j)                const;

     iterator                       column_secondary_end();
     const_iterator                 column_secondary_end()                 const;

     iterator                       column_secondary_end(size_type j);
     const_iterator                 column_secondary_end(size_type j)      const;

     reverse_iterator               column_secondary_rend();
     const_reverse_iterator         column_secondary_rend()                const;

     reverse_iterator               column_secondary_rend(size_type j);
     const_reverse_iterator         column_secondary_rend(size_type j)     const;

     reverse_iterator               column_rend();
     const_reverse_iterator         column_rend()                          const;

     reverse_iterator               column_rend(size_type j);
     const_reverse_iterator         column_rend(size_type j)               const;

     size_type                      size()                                 const;
     size_type                      linear_size()                          const;
     size_type                      rank()                                 const;
     size_type                      max_size()                             const;

     void                           resize(size_type new_size);
     void                           resize(size_type new_size, T c);

     size_type                      capacity()                             const;

     bool                           empty()                                const;

     void                           reserve(size_type n);

     reference                      operator() (size_type offset);
     const_reference                operator() (size_type offset)          const;

     reference                      operator() (size_type i, size_type j);
     const_reference                operator() (size_type i, size_type j)  const;

     reference                      operator[] (size_type offset);
     const_reference                operator[] (size_type offset)          const;

     reference                      at(size_type offset);
     const_reference                at(size_type offset)                   const;

     reference                      at(size_type i, size_type j);
     const_reference                at(size_type i, size_type j)           const;

     void                           clear();  
     void                           fill(const T& t);

     bool                           operator==(const TriangleMatrix<T,Allocator>& x);
     bool                           operator!=(const TriangleMatrix<T,Allocator>& x);
     bool                           operator< (const TriangleMatrix<T,Allocator>& x);
     bool                           operator>=(const TriangleMatrix<T,Allocator>& x);
     bool                           operator> (const TriangleMatrix<T,Allocator>& x);
     bool                           operator<=(const TriangleMatrix<T,Allocator>& x);

     // Added for use with lapack.
     //
           T*   raw_storage()       { return &my_vector[0]; }
     const T*   raw_storage() const { return &my_vector[0]; }

}; // end of class definition

#include "numerics/trimatrix.ipp"

} // end of namespace SAGE

#endif
