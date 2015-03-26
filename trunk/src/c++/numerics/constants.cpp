#include <iostream>
#include <iomanip>
#include <string>
#include "numerics/constants.h"

using namespace std;
using namespace SAGE;
typedef numeric_constants<double> constants;

#ifdef _MSC_VER

#undef max
#undef min

static const double _log_2 = log(2.0);

double log2(double x)
{
	return log(x)/_log_2;
}
#endif

void printval(ostream &out, const string& s, double n)
{
  out  << setw(15) << s << setw(30) << setprecision(24) << n;
  if(n>0)
    out << "\t" << "2^" << setprecision(5)  << log2(n);
  out << endl;
}

void printval2(ostream &out, const string& s, double n1, double n2)
{
  out  << setw(15) << s << setw(30) << setprecision(24) << n1
       << "\t" << setw(30) << setprecision(24)  << n2 << endl;
}

int main()
{
  printval(cout, "epsilon:",     constants::epsilon());
  printval(cout, "max:",         constants::max());
  printval(cout, "min:",         constants::min());
  printval(cout, "max exponent", constants::max_exponent);
  printval(cout, "min exponent", constants::min_exponent);
  printval(cout, "maxlog:",      constants::maxlog());
  printval(cout, "minlog:",      constants::minlog());
  printval(cout, "pi:",          constants::pi());
}
