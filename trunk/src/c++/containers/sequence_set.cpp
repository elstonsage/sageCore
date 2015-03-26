#include <iostream>
#include "sequence_set.h"
#include <cassert>

using namespace std;
using namespace SAGE;

template<class BAN, class IT>
void print(const BAN&b, IT begin, IT end)
{
  cout << b;
  for(; begin != end; ++begin)
    cout << *begin << " ";
  cout << endl;
}

int main()
{
  sequence_set<int> ss;


  ss.insert(6);
  ss.insert(5);
  ss.insert(4);
  ss.insert(3);
  ss.insert(2);
  ss.insert(1);

  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.sequence_begin(),  ss.sequence_end()  );
  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
  print("R-Sequential: ", ss.sequence_rbegin(), ss.sequence_rend() );

  ss.erase(5);
  ss.insert(9);
  ss.insert(0);
  ss.erase(9);
  ss.insert(5);
  ss.erase(0);
  ss.push_back(1);
  ss.insert(4);
  ss.pop_back();
  
  print("Ordered:      ", ss.ordered_begin(),   ss.ordered_end()   );
  print("Sequential:   ", ss.sequence_begin(),  ss.sequence_end()  );
  print("R-Ordered:    ", ss.ordered_rbegin(),  ss.ordered_rend()  );
  print("R-Sequential: ", ss.sequence_rbegin(), ss.sequence_rend() );

  assert( ss.sequence_find(1) == ss.find(1) );
  assert( ss.sequence_find(0) == ss.sequence_end() );
  assert( ss.ordered_end() == ss.find(0) );
}
