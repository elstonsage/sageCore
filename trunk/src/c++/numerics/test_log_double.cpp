#include <iostream>

using namespace std;

#include "numerics/log_double.h"

main()
{
  log_double a(1.), b(2.), c(3.), d(b);

  cout << a.get_double() << ' ' << b.get_double() << ' ' << c.get_double() << ' ' << d.get_double() << endl;

  log_double e(a+b), f(b*c);

  cout << e.get_double() << " X " << f.get_double() << endl;

  f *= d;

  cout << f.get_double() << endl;

  f /= b;

  cout << f.get_double() << endl;

  f += b;

  cout << f.get_double() << endl;

//  f += d * b - c;

  cout << f.get_double() << endl;

  // Test to see what the behaviors are for qnan

  double qnan = numeric_limits<double>::quiet_NaN();

  cout << SAGE::isnan(qnan) << ' ' << SAGE::isnan(log(qnan)) << ' ' << SAGE::isnan(exp(qnan))
       << endl;

  qnan = 1.0;

  cout << SAGE::isnan(qnan) << ' ' << SAGE::isnan(log(qnan)) << ' ' << SAGE::isnan(exp(qnan))
       << endl;

  cout << log_double(0.0).get_double() << endl;
  cout << (log_double(0.0) += 3.7).get_double() << endl;

}
