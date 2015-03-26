#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include "LSF/parse_ops.h"
#include "numerics/kahan.h"
#include "numerics/sinfo.h"
#include "numerics/corinfo.h"

using namespace std;
using namespace SAGE;


int main()
{
  cout << "sizeof(KahanAdder<double>) = " << sizeof(KahanAdder<double>) 
       << endl;

  size_t rank;

  cin >> rank >> ws;

  CorrelationInfo c(rank);

  while( !cin.eof() )
  {
    std::vector<double> data;

    for(size_t i = 0; i < rank; ++i)
    {
      string s;
      cin >> s;

      double d;

      if( !cin )
        d = std::numeric_limits<double>::quiet_NaN();
      else if( toUpper(s) == "INF" )
        d = std::numeric_limits<double>::infinity();
      else if( toUpper(s) == "-INF" )
        d = -std::numeric_limits<double>::infinity();
      else if( toUpper(s) == "QNAN" )
        d = std::numeric_limits<double>::quiet_NaN();
      else
        d = str2doub(s);

      data.push_back(d);
    }
    c += data;
  }

  cout << endl << "Covariance matrix: " << endl;

  for(size_t i = 0; i < rank; ++i)
  {
    for(size_t j = 0; j <= i; ++j)
      cout << setw(12) << c.covariance(i,j) << " ";
    cout << endl;
  }

  cout << endl << endl << "Correlation matrix: " << endl;

  for(size_t i = 0; i < rank; ++i)
  {
    for(size_t j = 0; j <= i ; ++j)
      cout << setw(12) << c.correlation(i,j) << " ";
    cout << endl;
  }
}
