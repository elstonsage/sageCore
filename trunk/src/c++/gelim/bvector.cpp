#ifdef __GNUG__
#define __PUT_STATIC_DATA_MEMBERS_HERE
#endif

#ifdef __KCC
#  include <vector>
#  define bit_vector vector<bool>
   using std::vector;
#else
#  include <bvector.h>
#endif

#include <assert.h>
#include <iostream.h>

bool subset(const bit_vector& lhs, const bit_vector& rhs)
{
  bit_vector::const_iterator i = lhs.begin();
  bit_vector::const_iterator j = rhs.begin();
  for( ; i != lhs.end() && j != rhs.end(); i++, j++)
    if((*i) && !(*j)) return 0;
  
  return 1;
}




#include "ftest.h"



main()
{
  bit_vector b(1000);

  for(int i = 0; i != 1000; i++)
    b[i] = rand() % 2;

  Test(b, 1000, 100);
  for(int i = 0; i < 10; ++i)
    cout << i << endl;

  vector<bool> v(10);

  v[1] |= true;

  cout << v[0] << v[1] << v[2] << endl;
}
