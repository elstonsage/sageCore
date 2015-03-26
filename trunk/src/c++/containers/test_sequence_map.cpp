#include <iostream>
#include "containers/sequence_map.h"
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
  sequence_map<int,int> ss;

  ss[6]=6;
  ss[5]=5;
  ss[4]=4;
  ss[3]=3;
  ss[2]=2;
  ss[1]=1;

  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.sequence_begin(),  ss.sequence_end()  );
//  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
//  print("R-Sequential: ", ss.sequence_rbegin(), ss.sequence_rend() );

  ss.erase(5);
  ss[9]=9;
  ss[0]=0;
  ss.erase(9);
  ss[5]=5;
  ss.erase(0);
  ss[1]=1;
  ss[4]=4;
  ss.pop_back();
  
  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.sequence_begin(),  ss.sequence_end()  );
//  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
//  print("R-Sequential: ", ss.sequence_rbegin(), ss.sequence_rend() );

  assert( ss.sequence_find(1) == ss.ordered_find(1) );
  assert( ss.sequence_find(0) == ss.sequence_end() );
  assert( ss.ordered_end() == ss.ordered_find(0) );
}
