#ifndef HASH_INHERITANCE_H
#define HASH_INHERITANCE_H

#include "containers/bitfield.h"

inline size_t hash_inheritance_xor(const bit_field& b)
{
  size_t t = 0;

  for(size_t i = 0; i < b.size(); i += 16)
  {
    t ^= b.to_u16(i);
  }
  
  return t;
}

inline size_t hash_inheritance_sum(const bit_field& b)
{
  size_t t = 0;

  for(size_t i = 0; i < b.size(); i += 16)
  {
    t += b.to_u16(i);
  }

  return t;
}

struct hash1
{
  size_t operator()(const bit_field& b) const { return hash_inheritance_xor(b); }
};

struct hash2
{
  size_t operator()(const bit_field& b) const { return hash_inheritance_sum(b); }
};

#endif
