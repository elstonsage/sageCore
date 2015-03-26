// =======================
// Inline Member Functions
// =======================

// ========================
// const_lvector_iterator
// ========================

inline const_lvector_iterator::const_reference
     const_lvector_iterator::operator*() const
{ return *p; }

inline const_lvector_iterator::const_reference
    const_lvector_iterator::operator[](difference_type i) const 
{ return *(*this + i); }

inline const_lvector_iterator& const_lvector_iterator::operator++()
{ ++p; return *this; }

inline const_lvector_iterator const_lvector_iterator::operator++(int)
{ const_iterator tmp = *this; ++p; return tmp; }

inline const_lvector_iterator& const_lvector_iterator::operator--()
{ --p; return *this; }

inline const_lvector_iterator const_lvector_iterator::operator--(int) 
{ const_iterator tmp = *this; --p; return tmp; }

inline const_lvector_iterator&
    const_lvector_iterator::operator+=(difference_type i)
{ p += i; return *this; }

inline const_lvector_iterator&
    const_lvector_iterator::operator-=(difference_type i)
{ *this += -i; return *this; }

inline const_lvector_iterator const_lvector_iterator::operator+(difference_type i) const
{ const_iterator tmp = *this; return tmp += i; }

inline const_lvector_iterator const_lvector_iterator::operator-(difference_type i) const
{ const_iterator tmp = *this; return tmp -= i; }

inline const_lvector_iterator::difference_type
    const_lvector_iterator::operator-(const_iterator x) const 
{ return p - x.p; }

inline bool const_lvector_iterator::operator==(const const_iterator& x) const 
{ return p == x.p; }

inline bool const_lvector_iterator::operator!=(const const_iterator& x) const
{ return p != x.p; }

inline bool const_lvector_iterator::operator<(const const_iterator& x) const
{ return p < x.p; }

// ==========================
// lvector Inline Functions
// ==========================

inline Likelihood_Vector::Likelihood_Vector()
    : storage(), first(0), f(false), bits(0), max_bits(0), my_scale(1.0), my_valid(true)
{ fixed = 0; }

inline Likelihood_Vector::Likelihood_Vector(size_type n)
    : first(0), f(false), bits(n), max_bits(n), my_scale(1.0), my_valid(true)
{ storage.resize(size(), 0); fixed = (1 << n) - 1; }

inline Likelihood_Vector::Likelihood_Vector(const lvector& l)
    : first(l.first), fixed(l.fixed), f(l.f), bits(l.bits), max_bits(l.bits), my_scale(l.my_scale)
{
  storage = l.storage;
  my_valid = l.is_valid();
}

inline Likelihood_Vector::~Likelihood_Vector() { deallocate(); }

inline Likelihood_Vector& Likelihood_Vector::operator= (const lvector& rhs)
{
  if(this == &rhs) return *this;

  if(!rhs.is_valid())
  {
    set_valid(false);
    return *this;
  }

  bits = rhs.bits;

  if(rhs.size() > capacity()) max_bits = rhs.bits;

  // This will allocate more memory if necessary.  No need to call initialize
  storage = rhs.storage;

  first    = rhs.first;
  fixed    = rhs.fixed;
  f        = rhs.f;
  my_scale = rhs.my_scale;
  
  set_valid(rhs.is_valid());

  return *this;
}

/* Removed until implementation available that makes sense.
inline void Likelihood_Vector::compress()
{
  if(bits == max_bits) return;
  
  storage.resize( size() );
}
*/

inline Likelihood_Vector& Likelihood_Vector::operator*= (log_double d)
{
  my_scale *= d;

  return *this;
}

inline void Likelihood_Vector::flatten(long l)
{
  first    = 0;
  f        = true;
  bits     = l;
  fixed    = 0;
  my_scale = log_double(1.0);

  set_valid(true);

  resize(l);

  double d = 1.0 / size();

  for(size_t i = 0; i < size(); ++i)
    storage[i] = d;
}

inline Likelihood_Vector::const_reference Likelihood_Vector::operator[] (size_type n) const
{ return *(begin() + n); }

inline Likelihood_Vector::const_iterator Likelihood_Vector::begin() const
{ return const_iterator(&*storage.begin()); }

inline Likelihood_Vector::const_iterator Likelihood_Vector::end  () const
{ return const_iterator(&*storage.begin()) + size(); }

inline Likelihood_Vector::const_reference Likelihood_Vector::front() const { return *begin();     }
inline Likelihood_Vector::const_reference Likelihood_Vector::back () const { return *(end() - 1); }

inline Likelihood_Vector::size_type Likelihood_Vector::size() const
{ return 1 << bits; }

inline Likelihood_Vector::size_type Likelihood_Vector::capacity() const 
{ return 1 << max_bits; }

inline Likelihood_Vector::size_type Likelihood_Vector::bit_count() const 
{ return bits; }

inline Likelihood_Vector::size_type Likelihood_Vector::bit_capacity() const 
{ return max_bits; }

inline log_double Likelihood_Vector::log_scale() const
{ return my_scale; }

inline log_double Likelihood_Vector::log_total() const
{ return my_scale * total(); }


inline void Likelihood_Vector::swap(lvector& rhs)
{ 
  storage.swap(rhs.storage);
  std::swap(storage,  rhs.storage);
  std::swap(first,    rhs.first);
  std::swap(fixed,    rhs.fixed);
  std::swap(f,        rhs.f);
  std::swap(bits,     rhs.bits);
  std::swap(max_bits, rhs.max_bits);
  std::swap(my_scale, rhs.my_scale);
}

inline bool Likelihood_Vector::operator== (const lvector& rhs) const
{
  if(&rhs == this) return true;
  if(size() != rhs.size()) return false;

  if(my_scale != rhs.my_scale) return false;

  return (::memcmp(&*storage.begin(), &*rhs.storage.begin(), size() * sizeof(likelihood)) == 0);
}

inline bool Likelihood_Vector::operator!= (const lvector& rhs) const
{ return !(*this == rhs); }

inline void Likelihood_Vector::increment_value(equivalence_class e, equivalence_class d, likelihood p)
{
  storage[e] += p;

  if(storage[e] != 0 && fixed)
  {
    if(f)
    {
      fixed &= ~(first ^ e);
    }
    else
    {
      f = true;
      first = e;
    }
  }
}

inline Likelihood_Vector::equivalence_class Likelihood_Vector::start() const
{ equivalence_class i = (equivalence_class) -1; return bump_up(i); }

inline Likelihood_Vector::equivalence_class&
    Likelihood_Vector::bump_up(equivalence_class& i) const
{ 
  ++i; 
  while( (i ^ first) & fixed ) i += (i ^ first) & fixed; 

  // kbj:  Fix to force i to at most size()
  if(i > size())
    i=size();

  return i;
}

inline void Likelihood_Vector::deallocate()
{
  storage.resize(0);
}

inline void Likelihood_Vector::clear_bits()
{
  for(size_type i = 0; i < storage.size(); ++i)
    storage[i] = 0.0;
}

inline void Likelihood_Vector::resize(size_type n)
{
  if(n > max_bits) max_bits = n;

  bits = n;
  
  storage.resize(size());
}

inline void Likelihood_Vector::resize_and_clear_bits(size_type n)
{
  if(n > max_bits)
  {
    clear_bits();

    bits = max_bits = n;
    
    storage.resize(size(), 0);
  }
  else
  {
    bits = n;

    storage.resize(size());
    
    clear_bits();
  }
}
inline bool Likelihood_Vector::is_valid() const
{
  return my_valid;
}
inline bool Likelihood_Vector::set_valid(bool v)
{
  return my_valid = v;
}

