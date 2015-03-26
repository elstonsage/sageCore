// =======================
// Inline Member Functions
// =======================

inline fft_bit_count::fft_bit_count() : built(false), mm(NULL) { }

inline fft_bit_count::~fft_bit_count() { deallocate(); }

inline bool fft_bit_count::build(long bits)
{
  if(bits <= 0) return false;
  
  initialize(bits);
  
  return built = true;
}

inline meiosis_map* fft_bit_count::get_meiosis_map() const
{ return mm; }

inline long fft_bit_count::operator[](equivalence_class e) const
{ if(built) return storage[e]; else return -1; }

inline fft_bit_count::equivalence_class fft_bit_count::capacity() const
{ return storage.capacity(); }

inline fft_bit_count::equivalence_class fft_bit_count::member_count() const
{ return storage.size(); }

inline void fft_bit_count::deallocate()
{
  storage.resize(0);
}

inline void fft_bit_count::initialize(equivalence_class n)
{
  n = 1 << n;
  storage.resize(n);
}

