#include <iostream>
#include <cassert>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <time.h>
#include "containers/bit_remap.h"

using namespace std;

int main()
{
  vector<int> br(32);

  for(int i=0; i < 32; ++i)
    br[i] = i;

  int bits = 4;
  int reps = 2000000;

  vector<double> v1(1<<bits);
  vector<double> v2(1<<bits);
  vector<double> o1(1<<bits);
  vector<double> o2(1<<bits);

//random_shuffle(br.begin(), br.end());

  bit_remap<unsigned int> remap(br.begin(), br.end());

  cout << "Building initial vectors" << endl;
  for(int i = 0; i < v1.size(); ++i)
  {
    v1[i] = (double)i/v1.size();
    v2[i] = (double)(v1.size()-i)/v1.size();
  }

  cout << "sizeof(remap) = " << sizeof(remap)      << endl
       << "value_size    = " << remap.value_size   << endl
       << "char_size     = " << remap.char_size    << endl
       << "vector_size   = " << remap.vector_size  << endl
       << "storage_size  = " << remap.storage_size << endl;

  cout << "Starting base loop" << endl;
  clock_t start = clock();
  unsigned int i;

  for(int r=0; r < reps; ++r)
    for(i = 0; i < v1.size(); ++i)
      o1[i] = v1[i]*v2[i];

  clock_t loop_overhead = clock() - start;

  cout << "Starting remapped loop" << endl;

  start = clock();

  for(int r=0; r < reps; ++r)
    for(i = 0; i < v1.size(); ++i)
      o2[i] = v1[i]*v2[remap(i)];

  clock_t remap_time = clock() - start;

  for(i = 0; i < v1.size(); ++i)
    assert( o1[i] == o2[i] );

  cout << " base loop time: " << loop_overhead/1e6 << endl;
  cout << "     remap time: " << remap_time/1e6                   << endl;
  cout << " remap overhead: " << (remap_time-loop_overhead)/1e6   << endl;
  cout << "       slowdown: " << (double)remap_time/loop_overhead << endl;
}
