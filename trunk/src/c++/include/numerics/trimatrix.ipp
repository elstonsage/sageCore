//----------------------------------------------------------------------------
// Implementation of iterator_base
//----------------------------------------------------------------------------

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator_base::iterator_base()
{
     my_row              = 0;
     my_column           = 0;            
	       
     my_row_direction    = 1;
     my_column_direction = 1;
	       
     my_row_major        = true;

     my_shape            = LOWER_TRIANGLE;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator_base::iterator_base(size_type r, size_type c)
{
     my_row              = r;
     my_column           = c;
	       
     my_row_direction    = 1;
     my_column_direction = 1;
	       
     my_row_major        = true;

     my_shape            = LOWER_TRIANGLE;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator_base::iterator_base(size_type r, size_type c, int r_d, int c_d, bool r_m)
{
     my_row              = r;
     my_column           = c;          
	       
     my_row_direction    = r_d;
     my_column_direction = c_d;
	       
     my_row_major        = r_m;

     my_shape            = LOWER_TRIANGLE;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator_base::iterator_base(size_type r, size_type c, int r_d, int c_d, bool r_m, matrix_shape shp)
{
     my_row              = r;
     my_column           = c;          
	       
     my_row_direction    = r_d;
     my_column_direction = c_d;
	       
     my_row_major        = r_m;

     my_shape            = shp;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator_base::~iterator_base()
{ }

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::size_type
TriangleMatrix<T,Allocator>::iterator_base::row() const
{
     return my_row;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::size_type
TriangleMatrix<T,Allocator>::iterator_base::column() const
{
     return my_column;
}

template <class T, class Allocator>
inline
int
TriangleMatrix<T,Allocator>::iterator_base::row_direction() const
{
     return my_row_direction;
}

template <class T, class Allocator>
inline
int
TriangleMatrix<T,Allocator>::iterator_base::column_direction() const
{
     return my_column_direction;
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T,Allocator>::iterator_base::row_major() const
{
     return my_row_major;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::matrix_shape
TriangleMatrix<T,Allocator>::iterator_base::shape() const
{
     return my_shape;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_row(size_type row)
{
     my_row = row;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_column(size_type column)
{
     my_column = column;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_row_direction(int row_d)
{
     my_row_direction = row_d;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_column_direction(int column_d)
{
     my_column_direction = column_d;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_row_major(bool row_major)
{
     my_row_major = row_major;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::set_shape(matrix_shape shape)
{
     my_shape = shape;
}
 
template <class T, class Allocator>
inline
void
TriangleMatrix<T,Allocator>::iterator_base::advance_iterator(size_type rank, bool forward)
{
     bool      swap_flag   = false;
     bool      r_major     = my_row_major;
     size_type r           = std::max(my_row, my_column);
     size_type c           = std::min(my_row, my_column);
     int       r_d         = my_row_direction;
     int       c_d         = my_column_direction;
     size_type over        = r + 1;     
                  
     if( !forward )
     {
       r_d *= -1;
       c_d *= -1;
     }
                  
     if( my_shape == FULL )
       over = rank;
       
     if(     my_shape == UPPER_TRIANGLE
         || (my_shape == FULL && r != my_row) )
     {
       swap_flag = true;
       std::swap( r_d, c_d );
       r_major = !r_major;
     }
     
     if( r_major )
     {
       c += c_d;
       
       // Underflow handle.
       if( c > rank + 1 )
       {
         r += r_d;
         
         if( my_shape != FULL )
           c = r;
         else               
           c = rank - 1;
       }
       // Overflow handle.
       else if( c == over )
       {
         r += r_d;
         c  = 0;
       }
     }
     else
     {
       if( !forward && r >= rank )
         std::swap( r, c );

       r += r_d;
       
       // Overflow handle.
       if( r == rank )
       {
         c += c_d;
         
         if( my_shape != FULL )
           r = c;
         else
           r = 0;
       }
       // Underflow handle.
       else if( r > rank || (r < c && my_shape != FULL) )
       {
           c += c_d;
           r  = rank - 1;
       }
     }
     
     if( swap_flag )
       std::swap( r, c );
      
     my_row    = r;
     my_column = c;
}

//----------------------------------------------------------------------------
// Implementation of iterator
//----------------------------------------------------------------------------
template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator::iterator()
{ 
     my_matrix = NULL;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator::iterator(TriangleMatrix<T,Allocator>* m, size_type r, size_type c) : iterator_base(r,c)
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator::iterator(TriangleMatrix<T,Allocator>* m, size_type r, size_type c,
                                                int r_d, int c_d, bool r_m)
                                      : iterator_base(r, c, r_d, c_d, r_m)
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator::iterator(TriangleMatrix<T,Allocator>* m, size_type r, size_type c,
                                                int r_d, int c_d, bool r_m,  matrix_shape shp)
                                      : iterator_base(r, c, r_d, c_d, r_m, shp)
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::iterator::~iterator()
{ }

template <class T, class Allocator>
inline
bool
TriangleMatrix<T,Allocator>::iterator::operator==(const iterator_base& x) const            
{
     return iterator_base::row() == x.row() && iterator_base::column() == x.column();      
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T,Allocator>::iterator::operator!=(const iterator_base& x) const
{
     return !(*this == x);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::iterator&
TriangleMatrix<T,Allocator>::iterator::operator++()
{
     advance_iterator( my_matrix->rank(), true );
     return *this;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::iterator&
TriangleMatrix<T,Allocator>::iterator::operator--()
{
     advance_iterator( my_matrix->rank(), false );
     return *this;
} 

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::iterator
TriangleMatrix<T,Allocator>::iterator::operator++(int)
{
     iterator temp(*this);
     advance_iterator( my_matrix->rank(), true );
     return temp;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::iterator
TriangleMatrix<T,Allocator>::iterator::operator--(int)
{
     iterator temp(*this);
     advance_iterator( my_matrix->rank(), false );
     return temp;
} 


template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::pointer
TriangleMatrix<T,Allocator>::iterator::operator->() const
{
     return &( operator*() );
}	       
 
template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::reference
TriangleMatrix<T,Allocator>::iterator::operator*() const
{
     return my_matrix(iterator_base::row(), iterator_base::column());
}

//----------------------------------------------------------------------------
// Implementation of const_iterator
//----------------------------------------------------------------------------
template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::const_iterator::const_iterator()
{ 
     my_matrix           = NULL;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::const_iterator::const_iterator(const TriangleMatrix<T,Allocator>* m, size_type r, size_type c) : iterator_base(r,c)
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::const_iterator::const_iterator(const TriangleMatrix<T,Allocator>* m, size_type r, size_type c,
                                                            int r_d, int c_d, bool r_m)
                                            : iterator_base(r, c, r_d, c_d, r_m)                                
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::const_iterator::const_iterator(const TriangleMatrix<T,Allocator>* m, size_type r, size_type c,
                                                            int r_d, int c_d, bool r_m,  matrix_shape shp)
                                            : iterator_base(r, c, r_d, c_d, r_m, shp)                                
{
     my_matrix = m;
}

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>::const_iterator::~const_iterator()
{ }

template <class T, class Allocator>
inline
bool
TriangleMatrix<T,Allocator>::const_iterator::operator==(const iterator_base& x) const
{
     return iterator_base::row() == x.row() && iterator_base::column() == x.column(); 
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T,Allocator>::const_iterator::operator!=(const iterator_base& x) const
{
     return !(*this == x);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_iterator&
TriangleMatrix<T,Allocator>::const_iterator::operator++()
{
     advance_iterator( my_matrix->rank(), true );
                       
     return *this;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_iterator&
TriangleMatrix<T,Allocator>::const_iterator::operator--()
{
     advance_iterator( my_matrix->rank(), false );
                       
     return *this;
} 

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_iterator
TriangleMatrix<T,Allocator>::const_iterator::operator++(int)
{
     const_iterator temp(*this);
     advance_iterator( my_matrix->rank(), true );
     return temp;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_iterator
TriangleMatrix<T,Allocator>::const_iterator::operator--(int)
{
     const_iterator temp(*this);
     advance_iterator( my_matrix->rank(), false );
     return temp;
} 

template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_pointer
TriangleMatrix<T,Allocator>::const_iterator::operator->() const
{
     return &( operator*() );
}	       
 
template <class T, class Allocator>
inline
typename TriangleMatrix<T,Allocator>::const_reference
TriangleMatrix<T,Allocator>::const_iterator::operator*() const
{
     return my_matrix(iterator_base::row(), iterator_base::column());
}

//----------------------------------------------------------------------------
// Implementation of TriangleMatrix
//----------------------------------------------------------------------------

template <class T, class Allocator>
inline
TriangleMatrix<T, Allocator>::TriangleMatrix(const Allocator& alloc) : my_vector(alloc)
{
     my_rank = 0;
}

template <class T, class Allocator>
inline
TriangleMatrix<T, Allocator>::TriangleMatrix(size_type n)
{
     resize(n);
}

template <class T, class Allocator>
inline
TriangleMatrix<T, Allocator>::TriangleMatrix(size_type n, const T& value,
                               const Allocator& alloc) : my_vector(alloc)
{
     resize(n, value);
}

template <class T, class Allocator>
inline
TriangleMatrix<T, Allocator>::TriangleMatrix(const TriangleMatrix<T,Allocator>& x)
{
     my_rank   = x.my_rank;
     my_vector = x.my_vector;
}

template <class T, class Allocator>
inline
TriangleMatrix<T, Allocator>::~TriangleMatrix()
{ }

template <class T, class Allocator>
inline
TriangleMatrix<T,Allocator>&
TriangleMatrix<T, Allocator>::operator=(const TriangleMatrix<T,Allocator>& x)
{
     my_rank   = x.my_rank;
     my_vector = x.my_vector;
    
     return *this;
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::assign(size_type n, const T& t)
{
     my_vector.assign(n, t);
}

template <class T, class Allocator>
template <class InputIterator>
inline
void
TriangleMatrix<T, Allocator>::assign(InputIterator first, InputIterator last)
{
     my_vector.assign(first, last);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_begin()
{
     return iterator(this, 0, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_begin() const
{
     return const_iterator(this, 0, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_begin(size_type i)
{
     return iterator(this, i, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_begin(size_type i) const
{
     return const_iterator(this, i, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_secondary_begin()
{
     return iterator(this, 0, 0, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_secondary_begin() const
{
     return const_iterator(this, 0, 0, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_secondary_begin(size_type i)
{
     return iterator(this, i, i, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_secondary_begin(size_type i) const
{
     return const_iterator(this, i, i, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rbegin()
{
     return iterator(this, rank() - 1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rbegin() const
{
     return const_iterator(this, rank() - 1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rbegin(size_type i)
{
     return iterator(this, i, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rbegin(size_type i) const
{
     return const_iterator(this, i, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_rbegin()
{
     return iterator(this, rank() - 1, rank() - 1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_rbegin() const
{
     return const_iterator(this, rank() - 1, rank() - 1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_rbegin(size_type i)
{
     return iterator(this, i, i, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_rbegin(size_type i) const
{
     return const_iterator(this, i, i, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_begin()
{
     return iterator(this, 0, 0, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_begin() const
{
     return const_iterator(this, 0, 0, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_begin(size_type j)
{
     return iterator(this, j, j, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_begin(size_type j) const
{
     return const_iterator(this, j, j, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_secondary_begin()
{
     return iterator(this, rank() - 1, rank() - 1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_secondary_begin() const
{
     return const_iterator(this, rank() - 1, rank() - 1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_secondary_begin(size_type j)
{
     return iterator(this, j, j, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_secondary_begin(size_type j) const
{
     return const_iterator(this, j, j, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rbegin()
{
     return iterator(this, rank() - 1, 0, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rbegin() const
{
     return const_iterator(this, rank() - 1, 0, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rbegin(size_type j)
{
     return iterator(this, rank() - 1, j, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rbegin(size_type j) const
{
     return const_iterator(this, rank() - 1, j, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_rbegin()
{
     return iterator(this, rank() - 1, rank() - 1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_rbegin() const
{
     return const_iterator(this, rank() - 1, rank() - 1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_rbegin(size_type j)
{
     return iterator(this, rank() - 1, j, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_rbegin(size_type j) const
{
     return const_iterator(this, rank() - 1, j, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_end()
{
     return iterator(this, rank(), 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_end() const
{
     return const_iterator(this, rank(), 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_end(size_type i)
{
     return iterator(this, i + 1, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_end(size_type i) const
{
     return const_iterator(this, i + 1, 0, 1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_secondary_end()
{
     return iterator(this, rank(), rank(), 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_secondary_end() const
{
     return const_iterator(this, rank(), rank(), 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::row_secondary_end(size_type i)
{
     return iterator(this, i + 1, i + 1, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::row_secondary_end(size_type i) const
{
     return const_iterator(this, i + 1, i + 1, 1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rend()
{
     return iterator(this, -1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rend() const
{
     return const_iterator(this, -1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rend(size_type i)
{
     return iterator(this, i - 1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_secondary_rend(size_type i) const
{
     return const_iterator(this, i - 1, 0, -1, 1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_rend()
{
     return iterator(this, -1, -1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_rend() const
{
     return const_iterator(this, -1, -1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::row_rend(size_type i)
{
     return iterator(this, i - 1, i - 1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::row_rend(size_type i) const
{
     return const_iterator(this, i - 1, i - 1, -1, -1, true);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_end()
{
     return iterator(this, rank(), rank(), 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_end() const
{
     return const_iterator(this, rank(), rank(), 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_end(size_type j)
{
     return iterator(this, j + 1, j + 1, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_end(size_type j) const
{
     return const_iterator(this, j + 1, j + 1, 1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_secondary_end()
{
     return iterator(this, -1, -1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_secondary_end() const
{
     return const_iterator(this, -1, -1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::iterator
TriangleMatrix<T, Allocator>::column_secondary_end(size_type j)
{
     return iterator(this, j - 1, j - 1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_iterator
TriangleMatrix<T, Allocator>::column_secondary_end(size_type j) const
{
     return const_iterator(this, j - 1, j - 1, 1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rend()
{
     return iterator(this, rank() - 1, rank(), -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rend() const
{
     return const_iterator(this, rank() - 1, rank(), -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rend(size_type j)
{
     return iterator(this, rank() - 1, j + 1, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_secondary_rend(size_type j) const
{
     return const_iterator(this, rank() - 1, j + 1, -1, 1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_rend()
{
     return iterator(this, rank() - 1, -1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_rend() const
{
     return const_iterator(this, rank() - 1, -1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reverse_iterator
TriangleMatrix<T, Allocator>::column_rend(size_type j)
{
     return iterator(this, rank() - 1, j - 1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reverse_iterator
TriangleMatrix<T, Allocator>::column_rend(size_type j) const
{
     return const_iterator(this, rank() - 1, j - 1, -1, -1, false);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::size_type
TriangleMatrix<T, Allocator>::linear_size() const
{
     return my_vector.size();
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::size_type
TriangleMatrix<T, Allocator>::size() const
{
     return my_rank;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::size_type
TriangleMatrix<T, Allocator>::rank() const
{
     return my_rank;
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::size_type
TriangleMatrix<T, Allocator>::max_size() const
{
     size_type c = my_vector.max_size();
     return floor( sqrt(2*c + 0.25) - 0.4999 );
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::resize(size_type n)
{
     my_rank = n;
     size_type number_of_elements = n * ( n + 1 ) / 2;
     my_vector.resize(number_of_elements);
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::resize(size_type n, T c)
{
     my_rank = n;
     size_type number_of_elements = n * ( n + 1 ) / 2;
     my_vector.resize(number_of_elements, c);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::size_type
TriangleMatrix<T, Allocator>::capacity() const
{
     size_type c = my_vector.capacity();
     return floor( sqrt(2*c + 0.25) - 0.4999 );
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::empty() const
{
     return my_vector.empty();
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::reserve(size_type n)
{
     size_type number_of_elements = n * ( n + 1 ) / 2;
     my_vector.reserve(number_of_elements);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reference
TriangleMatrix<T, Allocator>::operator()(size_type offset)
{
     return my_vector[offset];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reference
TriangleMatrix<T, Allocator>::operator()(size_type offset) const
{
     return my_vector[offset];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reference
TriangleMatrix<T, Allocator>::operator()(size_type i1, size_type i2)
{
     size_type i         = std::max(i1,i2);
     size_type j         = std::min(i1,i2);

     return my_vector[i * ( i + 1 ) / 2 + j];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reference
TriangleMatrix<T, Allocator>::operator()(size_type i1, size_type i2) const
{
     size_type i         = std::max(i1,i2);
     size_type j         = std::min(i1,i2);

     return my_vector[i * ( i + 1 ) / 2 + j];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reference
TriangleMatrix<T, Allocator>::operator[](size_type offset)
{
     return my_vector[offset];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reference
TriangleMatrix<T, Allocator>::operator[](size_type offset) const
{
     return my_vector[offset];
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reference
TriangleMatrix<T, Allocator>::at(size_type offset)
{
     return my_vector.at(offset);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reference
TriangleMatrix<T, Allocator>::at(size_type offset) const
{
     return my_vector.at(offset);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::reference
TriangleMatrix<T, Allocator>::at(size_type i, size_type j)
{
     size_type x         = std::max(i,j);
     size_type y         = std::min(i,j);

     return my_vector.at(x * ( x + 1 ) / 2 + y);
}

template <class T, class Allocator>
inline
typename TriangleMatrix<T, Allocator>::const_reference
TriangleMatrix<T, Allocator>::at(size_type i, size_type j) const
{
     size_type x         = std::max(i,j);
     size_type y         = std::min(i,j);

     return my_vector.at(x * ( x + 1 ) / 2 + y);
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::clear()
{
     my_vector.clear();
}

template <class T, class Allocator>
inline
void
TriangleMatrix<T, Allocator>::fill(const T& t)
{
     std::fill(my_vector.begin(), my_vector.end(), t);
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator==(const TriangleMatrix<T,Allocator>& x)
{
     return my_vector == x.my_vector;
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator!=(const TriangleMatrix<T,Allocator>& x)
{
     return !(this == x);
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator<(const TriangleMatrix<T,Allocator>& x)
{
     return my_vector < x.my_vector;
}
            
template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator>=(const TriangleMatrix<T,Allocator>& x)
{
     return !(this < x);
}

template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator>(const TriangleMatrix<T,Allocator>& x)
{
     return (x < this);
}
            
template <class T, class Allocator>
inline
bool
TriangleMatrix<T, Allocator>::operator<=(const TriangleMatrix<T,Allocator>& x)
{
     return !(x < this);
}

// end of TriangleMatrix implementation
