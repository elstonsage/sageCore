#ifndef UTIL_GET_MEM
#define UTIL_GET_MEM

#include <iostream>
#include <vector>
#include <string>

///
/// A utility for counting memory.
template<typename T>
struct MemCount
{
  static size_t go(const T & t) { return sizeof(T); }
};

///
/// A handy function.
template<typename T>
size_t get_mem(const T & t) { return MemCount<T>::go(t); }

/// Compiler macros to make it easier to write mem counters:
#define MEM_FRIEND(T) friend struct MemCount<T>;
#define MEM_COUNT_BEGIN(T) template<> struct MemCount<T> { static size_t go(const T & t)
#define MEM_COUNT_END };
#define MEM_DBG(x) std::cout << #x << " uses " << get_mem(x) << " bytes." << std::endl;

//==========================
// Specialized mem counters:
//==========================

// Strings:

MEM_COUNT_BEGIN(std::string)
{
  return t.size();
}
MEM_COUNT_END

// Vectors of doubles:
MEM_COUNT_BEGIN(std::vector<double>)
{
  return t.size();
}
MEM_COUNT_END

// Vectors of strings:
MEM_COUNT_BEGIN(std::vector<std::string>)
{
  size_t x = 0;
  for(size_t i = 0; i < t.size(); ++i)
    x += t[i].size();
  return x;
}
MEM_COUNT_END


#endif
