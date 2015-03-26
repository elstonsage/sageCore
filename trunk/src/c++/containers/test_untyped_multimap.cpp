#include "containers/UntypedMultimap.h"
#include <iostream>

class foo {};
class bar {};

int main()
{
  SAGE::UntypedMultimap m;

  m.insert(foo());
  m.insert(foo());
  m.insert(bar());
  m.insert(foo());

  for(SAGE::ObjConstIterator<foo> i = m.begin<foo>(); i != m.end<foo>(); ++i)
    std::cout << "found a foo!" << std::endl;

  for(SAGE::ObjConstIterator<bar> i = m.begin<bar>(); i != m.end<bar>(); ++i)
    std::cout << "found a bar!" << std::endl;


  return 0;
}
