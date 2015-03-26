#include "numerics/mt.h"
#include "numerics/functions.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace SAGE;

int main()
{
  cout << scientific << inv_chi_square_cdf( 1,0.001) << endl;
  cout << scientific << inv_chi_square_cdf( 1,1e-6)  << endl;
  cout << scientific << inv_chi_square_cdf(10,0.001) << endl;
  cout << scientific << inv_chi_square_cdf(10,1e-6)  << endl;

  MT mt(4357);

  cout << "Integer test: " << endl;

  for (int j=0; j<1000; ++j) 
  {
    cout << setw(10) << mt.uniform_integer() << " ";

    if (j%8==7)
      cout << endl;
  }
  cout << endl << endl;

  cout << "Real test: " << endl;
  mt.reseed(4357);
  for (int j=0; j<1000; ++j) 
  {
    cout << fixed << setw(10) << setprecision(8) << mt.uniform_real() << " ";

    if (j%8==7)
      cout << endl;
  }
  cout << endl;

}
