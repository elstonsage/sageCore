#include <iostream>
#include <functional>
#include "containers/muxvector.h"

using namespace std;
using namespace SAGE;

struct ModLess : public std::binary_function<int, int, bool>
{
  ModLess(int b = 2) : base(b) { }
  bool operator()(int i, int j) const { return (i % base) < (j % base); }
  int base;
};

int main()
{
  typedef    MultiplexVector<int, ModLess> MV;

  MV ms( ModLess(5) );

  for(size_t i = 0; i < 40; ++i)
    ms.insert(i);

  MV::type_iterator s;
  MV::value_iterator v;

  for(s = ms.type_begin(); s != ms.type_end(); ++s)
  {
    cout << s->first << ": ";
    for(v = s->second.begin(); v != s->second.end(); ++v)
      cout << *v << " ";
    cout << endl;
  }

  MV::iterator i;
  MV::iterator e = ms.end();

  int last = -1;

  for( i = ms.begin(); i != e; ++i )
  {
    if( i.type()->first != last )
    {
      last = i.type()->first;
      cout << endl << last << ": ";
    }

    cout << *i << " ";
  }
  cout << endl;
}

// Expected output:
// 0: 0 5 10 15 20 25 30 35
// 1: 1 6 11 16 21 26 31 36
// 2: 2 7 12 17 22 27 32 37
// 3: 3 8 13 18 23 28 33 38
// 4: 4 9 14 19 24 29 34 39
