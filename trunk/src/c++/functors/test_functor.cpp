#include "functors/functors.h"
#include <iostream>
#include <sstream>

void testIf()
{
  SAGE::FUNCTORS::If<> i;
  
  i.setIfFunctor(boost::lambda::constant(true));
  i.setThenFunctor(std::cout << boost::lambda::constant(1));
  i.setElseFunctor(std::cout << boost::lambda::constant(2));
  
  i();
  
  i.setIfFunctor(boost::lambda::constant(false));
  
  i();
}

void testFor()
{
  SAGE::FUNCTORS::For<> l;
  
  l.setInitFunctor     (boost::lambda::constant(0));
  l.setContinueFunctor (boost::lambda::_1 >= (size_t)10);
  l.setIterFunctor     (boost::lambda::_1 + 1);
  l.setBodyFunctor     (std::cout << boost::lambda::_1);

  l();
}

struct while_helper
{
  bool go() const { std::cout << i; --i; return true; }
  mutable int i;
};

void testWhile()
{
  SAGE::FUNCTORS::While<> w;

  while_helper h;
  h.i = 10;
  
  w.setContinueFunctor(boost::lambda::var(h.i) > 2);
  w.setBodyFunctor(boost::bind(&while_helper::go, &h));
    
  w();

  SAGE::FUNCTORS::DoWhile<> dw;
  
  h.i = 10;

  dw.setContinueFunctor(boost::lambda::var(h.i) > 2);
  dw.setBodyFunctor(boost::bind(&while_helper::go, &h));
    
  dw();
}

void testSequence()
{
  typedef SAGE::FUNCTORS::Sequence<int, boost::mpl::vector<int, int> > stype;
}

int main()
{
  testIf       ();
  testFor      ();
  testWhile    ();
  testSequence (); 

  std::cout << std::endl;

  return 0;
}

