#ifndef HASH_FUN_H
#define HASH_FUN_H

#include <limits>

namespace SAGE {

template <class Key> struct hash { };

inline size_t __def_hash_string(const char* s)
{
  unsigned long h = 0; 
  for ( ; *s; ++s)
    h = 5*h + *s;
  
  return size_t(h);
}

template<> struct hash<char*>
{
  size_t operator()(const char* s) const { return __def_hash_string(s); }
};

template<> struct hash<const char*>
{
  size_t operator()(const char* s) const { return __def_hash_string(s); }
};

template<> struct hash<char> {
  size_t operator()(char x) const { return x; }
};
template<> struct hash<unsigned char> {
  size_t operator()(unsigned char x) const { return x; }
};
template<> struct hash<signed char> {
  size_t operator()(unsigned char x) const { return x; }
};
template<> struct hash<short> {
  size_t operator()(short x) const { return x; }
};
template<> struct hash<unsigned short> {
  size_t operator()(unsigned short x) const { return x; }
};
template<> struct hash<int> {
  size_t operator()(int x) const { return x; }
};
template<> struct hash<unsigned int> {
  size_t operator()(unsigned int x) const { return x; }
};
template<> struct hash<long> {
  size_t operator()(long x) const { return x; }
};
template<> struct hash<unsigned long> {
  size_t operator()(unsigned long x) const { return x; }
};

/// Converts a double into a size_t by using the frexp() function.
/// This is a fast way to convert a double to a size_t that has good
/// properties for hashing.
///
/// The frexp() calculates a mantissa, exponent pair such that:
///
/// 0.5 < | mantissa | <= 1.0
/// and d = mantissa * 2 ^ exponent.
///
/// There is a unique pair for any double != 0.0
template<>
struct hash<double>
{
  size_t operator()(double d)
  {
    if(d == 0.0) return 0;
    
    int exponent;
    
    double mantissa = frexp(d, &exponent);
    
    return (size_t) (std::numeric_limits<size_t>::max() * (2.0 * fabs(mantissa) - 1.0));
  }
};

}

#endif
