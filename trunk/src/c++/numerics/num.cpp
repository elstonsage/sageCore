// This is a test

/* DESCRIPTION:
 *
 * This file contains a number of mathematical constants and
 * also some needed size parameters of the computer arithmetic.
 *
 * The default size parameters are as follows.
 *
 * For IEEE arithmetic:
 * MACHEP =  1.11022302462515654042E-16       2**-53
 * MAXLOG =  7.09782712893383996843E2         log(2**1024)
 * MINLOG = -7.08396418532264106224E2         log(2**-1022)
 * MAXNUM =  1.7976931348623158E308           2**1024

 * The global symbols for mathematical constants are
 * PI     =  3.14159265358979323846           pi
 * PIO2   =  1.57079632679489661923           pi/2
 * PIO4   =  7.85398163397448309616E-1        pi/4
 * SQRT2  =  1.41421356237309504880           sqrt(2)
 * SQRTH  =  7.07106781186547524401E-1        sqrt(2)/2
 * LOG2E  =  1.4426950408889634073599         1/log(2)
 * SQ2OPI =  7.9788456080286535587989E-1      sqrt( 2/pi )
 * LOGE2  =  6.93147180559945309417E-1        log(2)
 * LOGSQ2 =  3.46573590279972654709E-1        log(2)/2
 * THPIO4 =  2.35619449019234492885           3*pi/4
 * TWOOPI =  6.36619772367581343075535E-1     2/pi
 *
 * These lists are subject to change.
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <string>

namespace SAGE
{

template <class T>
class numeric_constants : public std::numeric_limits<T>
{
};


template<>
class numeric_constants<double> : public std::numeric_limits<double>
{
public:
  // Typically 2^-1022 also DBL_MIN -- fixme
  static double underflow_threshold() { return min(); }
  static double maxlog()              { return (max_exponent-1.0)*log_2();    }
  static double minlog()              { return (min_exponent-digits)*log_2(); }
  static double pi()                  { return 3.14159265358979323846;      }
  static double pi_over2()            { return 1.57079632679489661923;      }
  static double pi_over4()            { return 7.85398163397448309616E-1;   }
  static double log2_e()              { return 1.4426950408889634073599;    }
  static double sqrt_2overpi()        { return 7.9788456080286535587989E-1; }
  static double log_2()               { return 6.93147180559945309417E-1;   }
  static double log_sqrt2()           { return 3.46573590279972654709E-1;   }
  static double three_pi_over_4()     { return 2.35619449019234492885;      }
  static double two_over_pi()         { return 6.36619772367581343075535E-1;}
  static double negative_zero()       { return -0.0;                        }
};

}

using namespace std;
using namespace SAGE;
typedef numeric_constants<double> constants;

#ifdef _MSC_VER

#undef max
#undef min

static const double _log_2 = log(2.0);

static double log2(double x)
{
	return log(x)/_log_2;
}
#endif

void printval(ostream &out, const string& s, double n)
{
  out  << setw(15) << s << setw(30) << setprecision(24) << n;
  if(n>0)
    out << "\t" << "2^"     << setprecision(5)  << log2(n);
  out << endl;
}

void printval2(ostream &out, const string& s, double n1, double n2)
{
  out  << setw(15) << s << setw(30) << setprecision(24) << n1
       << "\t" << setw(30) << setprecision(24)  << n2 << endl;
}

int main()
{
  printval(cout, "epsilon:", constants::epsilon());
  printval(cout, "underflow:" , constants::underflow_threshold());
  printval(cout, "max:", constants::max());
  printval(cout, "min:", constants::min());
  printval(cout, "max exponent", constants::max_exponent);
  printval(cout, "min exponent", constants::min_exponent);
  printval(cout, "maxlog:"    , constants::maxlog());

  cout << log(pow(2.0,constants::max_exponent-1)) << endl;
  cout << constants::min_exponent-53 << endl;

  printval(cout, "minlog:"    , constants::minlog());

  cout << log(pow(2.0,constants::min_exponent-53)) << endl;
  printval(cout, "pi:"    , constants::pi());

  double pi = constants::pi();
  printval2(cout, "pi_over2:"  , constants::pi_over2(), pi/2.0);
  printval2(cout, "pi_over4:"  , constants::pi_over4(), pi/4.0);
  printval2(cout, "log2_e:"    , constants::log2_e(), log2(exp(1.0)));
  printval2(cout, "sqrt(2/p)i:", constants::sqrt_2overpi(), sqrt(2.0/pi));
  printval2(cout, "log(2):"   , constants::log_2(), log(2.0));
  printval2(cout, "log(2)/2:"  , constants::log_sqrt2(), log(2.0)/2.0);
  printval2(cout, "3pi/4:"     , constants::three_pi_over_4(), 3.0*pi/4.0);
  printval2(cout, "2/pi:"      , constants::two_over_pi(), 2.0/pi);
  printval2(cout, "-zero:"     , constants::negative_zero(), -0.0);
}

