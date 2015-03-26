#ifndef BIT_REMAP
#define BIT_REMAP

//****************************************************************************
//* File:      bit_remap.h                                                   *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Version 0.1                                                   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <functional>
#include <vector>
#include "globals/config.h"
#include "containers/bitfield.h"

template <class T, class S = unsigned char>
class bit_remap : std::unary_function<T,T>
{
public:

  typedef      T                value_type;
  typedef      S                storage_type;

  enum 
  { 
    value_size        = sizeof(value_type), 
    value_bits        = value_size*8,
    char_size         = sizeof(storage_type),
    char_bits         = char_size*8,
    vector_size       = value_size / char_size,
    storage_size      = 1 << (8*char_size)
  };
  
  value_type      bit_maps[vector_size][storage_size];

  bit_remap(unsigned int max_bits = value_bits)
  {
    clear();
    for(int bit = 0; bit < max_bits; ++bit)
      remap_bit(bit, bit);
  }
  // CHECK THIS!
  bit_remap(const bit_remap& r)
  {
    memcpy(bit_maps, r.bit_maps, &bit_maps[vector_size-1][storage_size-1]-bit_maps);
  }
  
  template <class InputIterator>
  bit_remap(InputIterator begin, InputIterator end)
  {
    clear();
    
    int bit;
    for(bit = 0; begin != end && bit < value_bits; ++begin, ++bit )
      remap_bit(bit, *begin);

    for( ; bit < value_bits; ++bit )
      remap_bit(bit, bit);
  }

  void clear()
  {
    for(int i=0; i < vector_size; ++i)
      for(int j = 0; j < storage_size; ++j)
        bit_maps[i][j] = 0;
  }

  void remap_bit(unsigned int bit, unsigned int new_bit)
  {
    if( bit >= value_bits || new_bit >= value_bits)
      return;

    int            map      = bit/char_bits;
    int            map_bit  = bit - 8*map; 
    storage_type   map_mask = 1U << map_bit;
    value_type     mask     = 1U << new_bit;

    for(unsigned int i = 0; i < storage_size; ++i)
      if(i & map_mask)
        bit_maps[map][i] |= mask;
  }

  inline
  value_type   operator()(value_type v) const
  {
    value_type  r = 0;
    value_type  s = storage_size-1;

    for(int i = 0; i < vector_size; ++i)
    {
      r  |= bit_maps[i][v & s];
      v >>= char_bits;
    }

    return r;
  }

};

template <>
class bit_remap<bit_field, unsigned char> : std::unary_function<bit_field&,bit_field>
{
public:

  typedef unsigned int                size_type;
  typedef bit_field                   value_type;
  typedef unsigned char               storage_type;
  typedef pair<size_type, value_type> value_pair;
  typedef vector<value_type>          value_vector;

  enum
  {
    char_size         = sizeof(storage_type),
    char_bits         = char_size*8,
    storage_size      = 1 << (8*char_size)
  };

  static const size_type bad_value = (size_type) -1;

  static const __U32  max_ss_value = storage_size-1;

  struct byte_mappings
  {
    size_type    offset;
    value_vector bits;

    byte_mappings() : offset(bad_value), bits(storage_size) { }
  };

  typedef vector<byte_mappings>         bit_field_mappings;

  bit_remap(unsigned int bits)
    : value_size(bits / 8 + ((bool) bits % 8)),
      value_bits(bits),
      vector_size(value_size)
  {
    bit_maps.resize(vector_size);
    
    clear();
    for(size_t bit = 0; bit < bits; ++bit)
      remap_bit(bit, bit);
  }
  
  void clear()
  {
    for(size_t i=0; i < vector_size; ++i)
    {
      bit_maps[i].offset = bad_value;
      for(size_t j = 0; j < storage_size; ++j)
        bit_maps[i].bits[j] = bit_field();
    }
  }

  void remap_bit(unsigned int bit, unsigned int new_bit)
  {
    if( bit >= value_bits || new_bit >= value_bits)
      return;

    int            map      = bit/char_bits; // byte of bit
    int            map_bit  = bit - 8*map;   // bit of byte
    storage_type   map_mask = 1U << map_bit; // mask of bit 

    unsigned int   offset     = new_bit >> 5;        // Number of 32 bit values 
    unsigned int   bit_offset = bit_maps[map].offset;
    unsigned int   length     = bit_maps[map].bits[0].size();

    if(offset < bit_offset)
    {
      size_type word_count = bit_offset - offset;

      if(bit_offset == bad_value) word_count = 1;

      for(unsigned int i = 0; i < storage_size; ++i)
      {
        bit_field b(length + (word_count << 5), false);

        for(size_t j = 0; j < word_count; ++j)
          b.set_chunk(bit_maps[map].bits[i].chunk(j<<5),(j+word_count) << 5);

        bit_maps[map].bits[i] = b;
      }
      bit_maps[map].offset = offset;
    }
    else if(new_bit > length + (bit_offset<<5))
    {
      for(unsigned int i = 0; i < storage_size; ++i)
      {
        bit_maps[map].bits[i].resize((offset - bit_offset + 1)<<5, false);
      }
      offset = bit_offset;
    }

    for(unsigned int i = 0; i < storage_size; ++i)
      if(i & map_mask)
      {
        bit_maps[map].bits[i][new_bit - (offset<<5)] = true;
      }
  }

  inline
  value_type   operator()(value_type v) const
  {
    value_type  r(value_bits, false);

    for(size_t i = 0; i < vector_size; )
    {
      __U32 u = v.to_u32(i * 8);
      
      size_t end = std::min(i+4, (size_t) vector_size);

      for( ; i < end; ++i)
      {
        size_type offset = bit_maps[i].offset;
        size_type length = bit_maps[i].bits[0].size();
        
        if(offset == bad_value) continue;

        for(size_type j = 0; j < length; j+=32)
          r.or_chunk(bit_maps[i].bits[u & max_ss_value].chunk(j),j+(offset<<5));

        u >>= char_bits;
      }
    }

    return r;
  }

private:

  bit_field_mappings bit_maps;

  size_t  value_size;
  size_t  value_bits;
  size_t  vector_size;
};

#endif
