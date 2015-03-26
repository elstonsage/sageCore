#ifndef SAMPLEINFO_H
#define SAMPLEINFO_H

#include <limits>
#include <algorithm>
#include <cmath>
#include "numerics/isnan.h"
#include "globals/config.h"
#include "numerics/kahan.h"

namespace SAGE
{

// Compute the following basic statistics for a given sample x1, x2, ...xN:
//   sum      - sum of x
//   sum1     - "      "
//   sum2     - sum of x^2
//   sum3     - sum of x^3
//   sum4     - sum of x^4
//   mean     - mean of x
//   variance - variance of x
//   skewness - skewness of x
//   kurtosis - kurtosis of x
//   count    - N
//   min      - lower bound of x
//   max      - upper bound of x
//
// A small sample adjustment can be given to correct for effects of small
// sample size.  This is currently only used in the computation of the
// variance as the term N/(N-a) is used to adjust the variance, where a is
// the small sample adjustment.

class SampleInfo
{
public:

  SampleInfo() { clear(); }
  explicit SampleInfo(size_t a) { clear(); my_sample_adjustment = a; }

  //typedef double             float_type;
  typedef KahanAdder<double> float_type;

  double sum()  const { return x1; }
  double sum1() const { return x1; }
  double sum2() const { return x2; }
  double sum3() const { return x3; }
  double sum4() const { return x4; }

  double mean() const
  {
    update();
    return my_mean;
  }

  double variance() const
  {
    update();
    return my_variance;
  }

  double standard_deviation() const
  {
    update();
    if( !finite(my_variance) )
      return my_variance;

    if( my_variance < 0 )
      return std::numeric_limits<double>::quiet_NaN();

    return sqrt(my_variance);
  }

  double standard_error() const
  {
    update();
    if( !my_count )
      return std::numeric_limits<double>::quiet_NaN();

    double dev = standard_deviation();

    if(!finite(dev))
      return dev;

    return dev/sqrt( (double)my_count );
  }

  double variation_coefficient() const
  {
    update();
    return standard_deviation()/mean();
  }

  double skewness() const
  {
    update();
    return my_skewness;
  }

  double kurtosis() const
  {
    update();
    return my_kurtosis;
  }

  size_t count() const
  {
    return my_count;
  }

  double max() const
  {
    return my_max;
  }

  double min() const
  {
    return my_min;
  }

  size_t sample_adjustment() const
  {
    return my_sample_adjustment;
  }

  void sample_adjustment(size_t a)
  {
    my_sample_adjustment = a;
  }

  void clear()
  {
    // Init the Kahan adders for x^n ( 0<n<5 )
    x1 = x2 = x3 = x4 = 0;

    // Reset counters for number of items
    my_count = my_last_count = 0;

    my_sample_adjustment = 0;

    // Init the return values
    my_mean = 0;

    const double temp = std::numeric_limits<double>::quiet_NaN();

    my_variance = my_skewness = my_kurtosis = temp;

    const double temp1 = std::numeric_limits<double>::infinity();

    my_max = -temp1;
    my_min =  temp1;
  }

  void add(double y)
  {
    // NaNs are not considered part of the sample
    if(SAGE::isnan(y))
      return;

    ++my_count;

#if 0
    // Handle positive and negative infinity
    if(!finite(y))
    {
      x1 = y;
      x2 = fabs(y);
      x3 = y;
      x4 = fabs(y);

      // Check to see if y is +/- infinity
      if( y > 0 )
        my_max = y;
      else
        my_min = y;

      return;
    }
#endif

    my_max = std::max(my_max, y);
    my_min = std::min(my_min, y);

    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y2*y2;

    x1 += y;
    x2 += y2;
    x3 += y3;
    x4 += y4;
  }

  SampleInfo& operator+=(double x) { add(x); return *this; }

  SampleInfo& operator+=(const SampleInfo& s)
  {
    my_count += s.my_count;

    my_sample_adjustment = std::max(my_sample_adjustment, s.my_sample_adjustment);

    my_max = std::max(my_max, s.my_max);
    my_min = std::min(my_min, s.my_min);

    x1 += s.x1;       // xn: sum of x^n
    x2 += s.x2;
    x3 += s.x3;
    x4 += s.x4;

    return *this;
  }

  SampleInfo operator+(const SampleInfo& s) const
  {
    SampleInfo temp(*this);
    return temp += s;
  }

private:

  void update() const
  {
    if(my_count == my_last_count)
      return;

    my_last_count = my_count;

    my_mean     = std::numeric_limits<double>::quiet_NaN();
    my_variance = std::numeric_limits<double>::quiet_NaN();
    my_skewness = std::numeric_limits<double>::quiet_NaN();
    my_kurtosis = std::numeric_limits<double>::quiet_NaN();

    if( !finite(x4) )
      return;

    double nx1 = x1 / my_count;     // Normalize sum(x^n)
    double nx2 = x2 / my_count;
    double nx3 = x3 / my_count;
    double nx4 = x4 / my_count;

    double m1 = nx1;               // Moments
    double m2 = nx2 - nx1*nx1;
    double m3 = nx3 - 3*nx2*nx1 + 2*nx1*nx1*nx1;
    double m4 = nx4 - 4*nx3*nx1 + 6*nx2*nx1*nx1 - 3*nx1*nx1*nx1*nx1;

    my_mean = m1;

    if( m2 < 10*std::numeric_limits<double>::epsilon() )
      return;

    if( my_count > my_sample_adjustment && my_count > 1 )
      my_variance = (double)my_count/(my_count-my_sample_adjustment)*m2;

    if( my_count > 2 )
      my_skewness = m3/pow(m2,3.0/2.0);

    if( my_count > 3 )
      my_kurtosis = m4/(m2*m2);
  }

  size_t my_count;
  size_t my_sample_adjustment;

  float_type x1;       // xn: sum of x^n
  float_type x2;
  float_type x3;
  float_type x4;

  double my_min;
  double my_max;

  // Mutable cache of computed values and when they were last computed
  mutable size_t my_last_count;
  mutable double my_mean;
  mutable double my_variance;
  mutable double my_skewness;
  mutable double my_kurtosis;
};

}

#endif
