#include <iostream>
#include "containers/indexed_map.h"
#include <cassert>

using namespace std;
using namespace SAGE;

template<class BAN, class IT>
void print(const BAN&b, IT begin, IT end)
{
  cout << b;
  for(; begin != end; ++begin)
    cout << begin->second << " ";
  cout << endl;
}

int main()
{
  indexed_map<string,int> ss;

  ss["f"]=6;
  ss["e"]=5;
  ss["d"]=4;
  ss["c"]=3;
  ss["b"]=2;
  ss["a"]=1;

  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.indexed_begin(),  ss.indexed_end()  );
//  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
//  print("R-Sequential: ", ss.indexed_rbegin(), ss.indexed_rend() );

  ss.erase("e");
  ss["i"]=9;
  ss["Z"]=0;

  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.indexed_begin(),  ss.indexed_end()  );

  ss.erase("i");
  ss["e"]=5;
  ss.erase("Z");


  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.indexed_begin(),  ss.indexed_end()  );

  ss["a"]=2;
  ss["d"]=3;
  ss.pop_back();
  
  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.indexed_begin(),  ss.indexed_end()  );
//  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
//  print("R-Sequential: ", ss.indexed_rbegin(), ss.indexed_rend() );

  assert( ss.indexed_find("a") == ss.ordered_find("a") );
  assert( ss.indexed_find("Z") == ss.indexed_end() );
  assert( ss.ordered_end() == ss.ordered_find("Z") );
}
