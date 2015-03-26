#include "containers/Reporter.h"
#include "containers/Sequence.h"
#include <iostream>
#include <sstream>

class Int2Int { public: int operator() (int i) const { return i; } };

class Int2String { public: std::string operator() (int i) const { std::ostringstream s; s << "final " << i; return s.str(); } };

class Void2Int { public: int operator() () const { return 3; } };

int main()
{
  Int2Int i2i;
  Int2String i2s;
  
  SAGE::Sequence<int, boost::mpl::vector<int> > seq;
  
/*  
  SAGE::Sequence<void> s2;
  
  NoResult noresult;
  
  s2.addFunctor("noresult", noresult);
  
  s2();
*/

  return 0;
}
