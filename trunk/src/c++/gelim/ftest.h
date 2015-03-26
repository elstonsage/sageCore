//
//  Bit storage class testing suite.
//
//  History: 0.1 gcw Initial Implementation               Aug    1997
//           1.0 gcw Wrote more detailed version, split   Mar 31 1998
//                   up the tests into suite of tests
//
// Copyright (c) 1998 R. C. Elston

#include <iostream.h>

// The tests are as follows:
// 1.  Constructor testing - Copy contruction
// 2.  Setting bits
// 3.  Equality tests
// 4.  Iteration test
// 5.  Random Access test
// 6.  Operator test - Takes two classes and Operator class.  For use with
//                     +, -, <<, etc.

template <class T>
int Test1(T& i, int sz, int count)
{
  unsigned long k = 0;
  
  for(int c = 0; c < count; c++)
  {
    T a = i;
    a[0] = 1;  // To make sure a is constructed
    k++;
  }
  
  return k;
}

template <class T>
int Test2(T& i, int sz, int count)
{
  unsigned long k = 0;
  
  for(int c = 0; c < count; c++)
  {
    i[count % size] = !i[count % size];
    k++;
  }
  
  return k;
}

template <class T>
int Test3(T& i, int sz, int count)
{
  unsigned long k = 0;
  
  T a = i;
  bool b;
  for(int c = 0; c < count; c++)
  {
    i[count % size] = !i[count % size];
    b = (a == i);
    k++;
  }
  
  return k;
}

template <class T>
int Test4(T& i, int sz, int count)
{
  unsigned long k = 0;
  
  for(int c = 0; c < count; c++)
  {
    for(T::iterator it= i.begin(); it != i.end(); it++, k++);
  }
  
  return k;
}

template <class T>
int Test5(T& i, int sz, int count)
{
  unsigned long k = 0;
  
  for(int c = 0; c < count; c++)
  {
    i[rand() % size] = !i[rand() % size];
    k++;
  }
  
  return k;
}

template <class T, class S, class CMP>
int Test6(T& i, S& j, int sz, int count, CMP& cmp)
{
  unsigned long k = 0;
  
  for(int c = 0; c < count; c++)
  {
    cmp(i, j);
    k++;
  }
  
  return k;
}
