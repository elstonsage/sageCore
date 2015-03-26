#include "util/AutoTrace.h"

#include <boost/mpl/for_each.hpp>

class C
{
public:
  C();
  
  void doNothing() { }

  static void foo();
};

void
AUTOTRACE_STATIC(C, NOARGS,
C::foo() {)
}

class B
{
public:
  B();
  void doNothing() { }
};

class A
{
public:
  A(int a, int b);

  int go(int y);
private:

  int my_x1;
  int my_x2;
};

AUTOTRACE(NOARGS,
C::C() 
{)
}

AUTOTRACE(NOARGS,
B::B() 
{)
  C c; 
  
  c.doNothing();
}

int 
AUTOTRACE(ARGLIST1(y),
A::go(int y)
{)
  B b;
  
  b.doNothing();

  return 0;
}

AUTOTRACE_CT2(NOARGS, 
A::A(int a, int b) : my_x1(0), my_x2(1) 
{)

}

int main(int argc, char* argv[])
{
  SAGE::UTIL::AutoTrace::excludeClass<B> ();

  A a(3,4);

  a.go(4);

  C::foo();

  return 0;
}
