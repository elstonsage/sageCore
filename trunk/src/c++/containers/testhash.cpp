#include <stdlib.h>
#include <iostream>
#include "containers/cache_map.h"

using namespace std;
using namespace SAGE;

int main()
{
  typedef cache_map<size_t,size_t> map_type;

  srand( time(NULL) );
  
  map_type m;
  m.resize(4516100);


  cout << "Begin" << endl;

  for(size_t i = 0; i < 10000000; ++i)
  {
    int r = rand() % 10000000;
    m[r] = r;
  }

  map_type::const_iterator j = m.begin();
  
  for( ; j != m.end(); ++j)
    if(j->is_set() && (*j)->first != (*j)->second)
      cout << "Error! Mismatch1: " << (*j)->first  << " " 
                                   << (*j)->second << endl;

  for(size_t i = 0; i < 10000000; ++i)
  {
    map_type::const_iterator j = m.find(i);
    if(j != m.end() && (*j)->first != (*j)->second && (*j)->first == i )
      cout << "Error! Mismatch2: " << i            << " " 
                                   << (*j)->first  << " " 
                                   << (*j)->second << endl;
 
  }

  for(size_t i = 0; i < 50000000; ++i)
  {
    int k = rand();
    map_type::const_iterator j = m.find(k);
    if(j != m.end() && (*j)->first != (*j)->second && (*j)->first == k )
      cout << "Error! Mismatch2: " << k            << " " 
                                   << (*j)->first  << " " 
                                   << (*j)->second << endl;
  }

  cout << "End" << endl;
}
