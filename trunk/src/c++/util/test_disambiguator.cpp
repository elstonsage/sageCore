#include "util/disambiguator.h"
#include <cassert>
#include <iostream>
#include <iomanip>

DISAMBIGUATE(TYPE1, int);
DISAMBIGUATE(TYPE2, int);

struct foo
{
  void bar() const { ; }
};


DISAMBIGUATE(TYPE3, const foo);
DISAMBIGUATE(TYPE4, const foo);

void test_func(TYPE1 t1, TYPE2 t2)
{
  std::cout << t1 << " + " << t2 << " = " << (t1 + t2) << std::endl;
  
  t1() = 5;
}  

void test_func2(TYPE3 t3, TYPE4 t4)
{
  t3().bar();
  t4().bar();
}  

struct test_in_class
{
  DISAMBIGUATE(TYPE1, const int);
  DISAMBIGUATE(TYPE2, const int);

  test_in_class(TYPE1 t1, TYPE2 t2) { }
};

int main()
{
  int i1 = 3, i2 = 4;
  
  test_func(TYPE1(i1), TYPE2(i2));
  
  assert(i1 == 5);

  test_func2(TYPE3(foo()), TYPE4(foo()));
  
  test_in_class::TYPE1 t1(0);
  test_in_class::TYPE2 t2(1);
  
  test_in_class(t1, t2);

  return 0;
}
