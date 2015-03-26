#include <iostream>
#include <functional>
#include "containers/muxset.h"

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

  typedef MultiplexSet<int, ModLess, greater<int> > MS;

  MS ms( ModLess(5) );

  for( size_t i = 0; i < 40; ++i )
    ms.insert(i);

  MS::type_iterator t;
  MS::value_iterator v;

  for( t = ms.type_begin(); t != ms.type_end(); ++t )
  {
    cout << t->first << ": ";
    for( v = t->second.begin(); v != t->second.end(); ++v )
      cout << *v << " ";
    cout << endl;
  }
  cout << endl;

  cout << "ms.type_empty()    = " << ms.type_empty() << endl;
  cout << "ms.type_size()     = " << ms.type_size() << endl;
  cout << "ms.type_max_size() = " << ms.type_max_size() << endl;
  cout << "ms.type_count(3)   = " << ms.type_count(3) << endl;  
  cout << endl;
  cout << "ms.value_empty(3)    = " << ms.value_empty(3) << endl;  
  cout << "ms.value_size(3)     = " << ms.value_size(3) << endl;  
  cout << "ms.value_max_size(3) = " << ms.value_max_size(3) << endl;  
  cout << "ms.value_count(3)    = " << ms.value_count(3) << endl;  
  cout << endl;

  MS::iterator k = ms.lower_bound(37);
  cout << "lower_bound(37) = " << *k << endl;
  k = ms.upper_bound(37);
  cout << "upper_bound(37) = " << *k << endl;
  cout << endl;

  std::pair<MS::iterator, MS::iterator> ms_i_pair = ms.equal_range(33);
  cout << "equal_range(33) start = " << *(ms_i_pair.first) << endl;
  cout << "equal_range(33) end   = " << *(ms_i_pair.second) << endl;
  cout << endl;

  MS::iterator i;
  MS::iterator j;
  i = ms.begin();
  j = ms.find(7);
  cout << "*i = " << *i << endl;
  cout << "*j = " << *j << endl;
  cout << endl;

  cout << "ms.erase_value(42) = "<< ms.erase_value(42) << endl;
  ms.erase_value(i);
  ms.erase_value(i, j);

  i = ms.begin();
  j = ms.find(8);
  cout << "ms.erase_type(41)  = "<< ms.erase_type(41) << endl;
  ms.erase_type(i);
  ms.erase_type(j);
  i = ms.begin();
  j = ms.find(7);
  ms.erase_type(i, j);
  cout << endl;

  for( t = ms.type_begin(); t != ms.type_end(); ++t )
  {
    cout << t->first << ": ";
    for( v = t->second.begin(); v != t->second.end(); ++v )
      cout << *v << " ";
    cout << endl;
  }
  cout << endl;

  for( i = ms.begin(); i != ms.end(); ++i )
    cout << *i << " ";
  cout << endl;

  for( i = ms.end(); i != ms.begin(); )
  {
    --i;
    cout << *i << " ";
  }
  cout << endl;

}

// Expected output:
// 0: 35 30 25 20 15 10 5 0
// 1: 36 31 26 21 16 11 6 1
// 2: 37 32 27 22 17 12 7 2
// 3: 38 33 28 23 18 13 8 3
// 4: 39 34 29 24 19 14 9 4
