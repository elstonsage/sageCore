#include <iostream>
#include <functional>
#include "containers/muxmap.h"

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
  free(malloc(1));

  typedef MultiplexMap<int, int, ModLess, greater<int> > MM;

  MM mm( ModLess(5) );

  for( size_t i = 0; i < 40; ++i )
    mm.insert(make_pair(i%5, i));


  MM::type_iterator t;
  MM::value_iterator v;

  for( t = mm.type_begin(); t != mm.type_end(); ++t )
  {
    cout << t->first << ": ";
    for( v = t->second.begin(); v != t->second.end(); ++v )
      cout << *v << " ";
    cout << endl;
  }
  cout << endl;

  cout << "mm.type_empty()    = " << mm.type_empty() << endl;
  cout << "mm.type_size()     = " << mm.type_size() << endl;
  cout << "mm.type_max_size() = " << mm.type_max_size() << endl;
  cout << "mm.type_count(3)   = " << mm.type_count(3) << endl;  
  cout << endl;
  cout << "mm.value_empty(3)    = " << mm.value_empty(3) << endl;  
  cout << "mm.value_size(3)     = " << mm.value_size(3) << endl;  
  cout << "mm.value_max_size(3) = " << mm.value_max_size(3) << endl;  
  cout << "mm.value_count(3)    = " << mm.value_count(make_pair(3%5, 3)) << endl;  
  cout << endl;


  MM::iterator k = mm.lower_bound(37);
  cout << "lower_bound(37) = " << *k << endl;
  k = mm.upper_bound(37);
  cout << "upper_bound(37) = " << *k << endl;
  cout << endl;

  std::pair<MM::iterator, MM::iterator> mm_i_pair = mm.equal_range(33);
  cout << "equal_range(33) start = " << *(mm_i_pair.first) << endl;
  cout << "equal_range(33) end   = " << *(mm_i_pair.second) << endl;
  cout << endl;


  MM::iterator i;
  MM::iterator j;
  i = mm.begin();
  j = mm.find(make_pair(7%5, 7));
  cout << "*i = " << *i << endl;
  cout << "*j = " << *j << endl;
  cout << endl;

  cout << "mm.erase_value(42) = " << mm.erase_value(make_pair(42%5, 42)) << endl;
  mm.erase_value(i);

  i = mm.begin();
  j = mm.find(make_pair(8%5, 8));
  cout << "mm.erase_type(41)  = " << mm.erase_type(41) << endl;
  mm.erase_type(i);
  mm.erase_type(j);
  i = mm.begin();
  j = mm.find(make_pair(7%5, 7));
  mm.erase_type(i, j);
  cout << endl;

  for( t = mm.type_begin(); t != mm.type_end(); ++t )
  {
    cout << t->first << ": ";
    for( v = t->second.begin(); v != t->second.end(); ++v )
      cout << *v << " ";
    cout << endl;
  }
  cout << endl;

  for( i = mm.begin(); i != mm.end(); ++i )
    cout << *i << " ";
  cout << endl;

  for( i = mm.end(); i != mm.begin(); )
  {
    --i;
    cout << *i << " ";
  }
  cout << endl << endl;

  cout << "operator[2] = ";
  std::set<int, greater<int> > s = mm[2];
  std::set<int, greater<int> >::iterator s_i;
  for( s_i = s.begin(); s_i != s.end(); ++s_i )
    cout << *s_i << " ";
  cout << endl << endl;
}

// Expected output:
// 0: 35 30 25 20 15 10 5 0
// 1: 36 31 26 21 16 11 6 1
// 2: 37 32 27 22 17 12 7 2
// 3: 38 33 28 23 18 13 8 3
// 4: 39 34 29 24 19 14 9 4
