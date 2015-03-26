namespace SAGE {
namespace SAMPLING {

//===============================================================================================================
//===============================================================================================================
//===============================================================================================================
//
//			NON-CONST ITERATOR
//
//===============================================================================================================
//===============================================================================================================
//===============================================================================================================

//===========================================================================
// CONSTRUCTOR - Protected
//===========================================================================
inline 
FieldIterator::FieldIterator(
	FieldVectorType          * _masterlist, 
	FieldGroupType::iterator   _iter) 
	:
	masterlist (_masterlist),
	iter       (_iter)
{}

//===========================================================================
// COPY CONSTRUCTOR
//===========================================================================
inline
FieldIterator::FieldIterator(const FieldIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// operator* ()
//===========================================================================
inline FieldIterator::reference
FieldIterator::operator*()
{ 
  return (*masterlist)[iter->id];
}

//===========================================================================
// operator[] (...)
//===========================================================================
inline FieldIterator::reference
FieldIterator::operator[](difference_type i)
{ 
  return *(*this + i); 
}
        
//===========================================================================
// operator-> ()
//===========================================================================
inline FieldIterator::pointer 
FieldIterator::operator->()
{
  return &(*masterlist)[iter->id];
}

//===========================================================================
// operator++ () POST INCREMENT
//===========================================================================
inline FieldIterator::iterator
FieldIterator::operator++()
{ 
  increment();
  return *this;
}
  
//===========================================================================
// operator++ (...) PRE INCREMENT
//===========================================================================
inline FieldIterator::iterator
FieldIterator::operator++(int)
{ 
  iterator tmp = *this; 
  increment(); 
  return tmp; 
}

//===========================================================================
// operator-- () POST DECREMENT
//===========================================================================
inline FieldIterator::iterator
FieldIterator::operator--()
{ 
  decrement();
  return *this;
}

//===========================================================================
// operator-- (...) PRE DECREMENT
//===========================================================================
inline FieldIterator::iterator
FieldIterator::operator--(int)
{ 
  iterator tmp = *this; 
  decrement(); 
  return tmp; 
}

//===========================================================================
// operator- (...)
//===========================================================================
inline FieldIterator::difference_type 
FieldIterator::operator-(const iterator & other) const
{
  return iter - other.iter;
}

//===========================================================================
// operator== (...)
//===========================================================================
inline bool
FieldIterator::operator== (const iterator & other) const
{
  return iter == other.iter;
}

//===========================================================================
// operator!= (...)
//===========================================================================
inline bool 
FieldIterator::operator!= (const iterator & other) const
{
  return !(*this == other);
}

//===========================================================================
// operator< (...)
//===========================================================================
inline bool 
FieldIterator::operator< (const iterator & other) const
{
  return iter < other.iter;
}

//===========================================================================
// increment()
//===========================================================================
inline void 
FieldIterator::increment()
{
  iter++;
}

//===========================================================================
// decrement()
//===========================================================================
inline void 
FieldIterator::decrement()
{
  iter--;
}

//===============================================================================================================
//===============================================================================================================
//===============================================================================================================
//
//			CONST ITERATOR
//
//===============================================================================================================
//===============================================================================================================
//===============================================================================================================

//===========================================================================
// CONSTRUCTOR - Protected
//===========================================================================
inline 
FieldConstIterator::FieldConstIterator(
	const FieldVectorType          * _masterlist, 
	FieldGroupType::const_iterator   _iter) 
	:
	masterlist (_masterlist),
	iter       (_iter)
{}

//===========================================================================
// COPY CONSTRUCTOR #1
//===========================================================================
inline
FieldConstIterator::FieldConstIterator(const FieldConstIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// COPY CONSTRUCTOR #2
//===========================================================================
inline
FieldConstIterator::FieldConstIterator(const FieldIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// operator* ()
//===========================================================================
inline FieldConstIterator::const_reference
FieldConstIterator::operator*() const
{ 
  return (*masterlist)[iter->id];
}

//===========================================================================
// operator[] (...)
//===========================================================================
inline FieldConstIterator::const_reference
FieldConstIterator::operator[](difference_type i) const
{ 
  return *(*this + i); 
}
        
//===========================================================================
// operator-> ()
//===========================================================================
inline FieldConstIterator::const_pointer 
FieldConstIterator::operator->()
{
  return &(*masterlist)[iter->id];
}

//===========================================================================
// operator++ () POST INCREMENT
//===========================================================================
inline FieldConstIterator::const_iterator
FieldConstIterator::operator++()
{ 
  increment();
  return *this;
}
  
//===========================================================================
// operator++ (...) PRE INCREMENT
//===========================================================================
inline FieldConstIterator::const_iterator
FieldConstIterator::operator++(int)
{ 
  const_iterator tmp = *this; 
  increment(); 
  return tmp; 
}

//===========================================================================
// operator-- () POST DECREMENT
//===========================================================================
inline FieldConstIterator::const_iterator
FieldConstIterator::operator--()
{ 
  decrement();
  return *this;
}

//===========================================================================
// operator-- (...) PRE DECREMENT
//===========================================================================
inline FieldConstIterator::const_iterator
FieldConstIterator::operator--(int)
{ 
  const_iterator tmp = *this; 
  decrement(); 
  return tmp; 
}

//===========================================================================
// operator- (...)
//===========================================================================
inline FieldConstIterator::difference_type 
FieldConstIterator::operator-(const const_iterator & other) const
{
  return iter - other.iter;
}

//===========================================================================
// operator== (...)
//===========================================================================
inline bool
FieldConstIterator::operator== (const const_iterator & other) const
{
  return iter == other.iter;
}

//===========================================================================
// operator!= (...)
//===========================================================================
inline bool 
FieldConstIterator::operator!= (const const_iterator & other) const
{
  return !(*this == other);
}

//===========================================================================
// operator< (...)
//===========================================================================
inline bool 
FieldConstIterator::operator< (const const_iterator & other) const
{
  return iter < other.iter;
}

//===========================================================================
// increment()
//===========================================================================
inline void 
FieldConstIterator::increment()
{
  iter++;
}

//===========================================================================
// decrement()
//===========================================================================
inline void 
FieldConstIterator::decrement()
{
  iter--;
}

} // End namespace SAMPLING
} // End namespace SAGE
