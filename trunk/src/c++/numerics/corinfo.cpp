#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include "numerics/corinfo.h"

namespace SAGE {

void
SimpleCorrelationInfo::add(double a, double b)
{
    my_info1 += a;
    my_info2 += b;

    if( SAGE::isnan(a) || SAGE::isnan(b) )
      return;

    ++my_count;

    my_sums.add( a, b );
}

SimpleCorrelationInfo&
SimpleCorrelationInfo::operator+=(const SimpleCorrelationInfo &c)
{
    my_count += c.my_count;

    my_info1 += c.my_info1;
    my_info2 += c.my_info2;
    my_sums  += c.my_sums;

    return *this;
}

void
SimpleCorrelationInfo::do_update() const
{
    my_last_count = my_count;

    // Initial return values.
    // my_correlations.resize( my_sums.size() );

    correlation_type&                   corr = my_correlation;

    // Init results to default values
    corr.covariance  = std::numeric_limits<double>::quiet_NaN();
    corr.correlation = std::numeric_limits<double>::quiet_NaN();

        // Cross product will be NaN or infinite if either x or y are.
    if( !finite(my_sums.xy) || my_count == 0 )
      return;

    const double x_mean = my_sums.x / my_count;
    const double y_mean = my_sums.y / my_count;

    corr.covariance = my_sums.xy / my_count - x_mean * y_mean;

    if( my_count < 3 )
      return;

    const double x_variance  = my_sums.xx / my_count - x_mean * x_mean;
    const double y_variance  = my_sums.yy / my_count - y_mean * y_mean;
    const double xy_variance = x_variance * y_variance;

    if( xy_variance >= std::numeric_limits<double>::epsilon() )
      corr.correlation = corr.covariance / sqrt(xy_variance);

} // end of update().

//---------------------------------------------------------------------------------


void
CorrelationInfo::resize(size_t variates)
{
    // Reset counters for number of items.
    my_count = my_last_count = 0;

    my_variates.resize(0);
    my_variates.resize(variates);

    my_sums.resize(variates);
    my_correlations.resize(variates);
}

void
CorrelationInfo::add(const std::vector<double>& data, double weight)
{
    // Can't add with different variates count.
    if( size() != data.size() )
      return;

    ++my_count;

    // loop index i, j == data[i][j].
    size_t i = 0, j = 0;

    // Add variates to SampleInfo.
    for( i = 0; i < size(); ++i )
       my_variates[i] += data[i];

    // Add to my_correlation_sums.
    for( i = 0; i < size(); ++i )
    {
      if( SAGE::isnan(data[i]) )
        continue;

      for( j = 0; j <= i; ++j )
      {
        if( SAGE::isnan(data[j]) )
          continue;

        my_sums(i, j).add( data[i], data[j], weight );
      }
    }
}

void
CorrelationInfo::add(const std::vector<double>& data, const TriangleMatrix<double>& weight)
{
    // Can't add with different variates count.
    if( size() != data.size() )
      return;

    ++my_count;

    // loop index i, j == data[i][j].
    size_t i = 0, j = 0;

    // Add variates to SampleInfo.
    for( i = 0; i < size(); ++i )
       my_variates[i] += data[i];

    // Add to my_correlation_sums.
    for( i = 0; i < size(); ++i )
    {
      if( SAGE::isnan(data[i]) )
        continue;

      for( j = 0; j <= i; ++j )
      {
        if( SAGE::isnan(data[j]) )
          continue;

        my_sums(i, j).add( data[i], data[j], weight(i, j) );
      }
    }
}

void
CorrelationInfo::add(size_t i, double di, size_t j, double dj, double weight)
{
    // Add variates to SampleInfo.
    my_variates[i] += di;
    my_variates[j] += dj;

    if( SAGE::isnan(di) || SAGE::isnan(dj) )
      return;

    ++my_count;

    my_sums(i, j).add( di, dj, weight );
}


CorrelationInfo&
CorrelationInfo::operator+=(const CorrelationInfo &c)
{
    // Can't add with different variates count.
    if( size() != c.size() )
      return *this;

    my_count += c.my_count;

    // Add c.my_variates to my_variates.
    for( size_t i = 0; i < size(); ++i )
      my_variates[i] += c.my_variates[i];

    // Add c.my_sums to my_sums.
    for( size_t i = 0; i < size(); ++i )
      for( size_t j = 0; j <= i; ++j )
        my_sums(i, j) += c.my_sums(i, j);

    return *this;
}

void
CorrelationInfo::do_update() const
{
    my_last_count = my_count;

    // Initial return values.
    // my_correlations.resize( my_sums.size() );

    // Update my_correlations values.
    for(size_t i = 0; i < size(); ++i )
    {
      for(size_t j = 0; j <= i; ++j )
      {
        // Create local aliases for long names to improve readibility
        const correlation_computation_type&  sum = my_sums(i, j);
        correlation_type&                   corr = my_correlations(i, j);

        // Init results to default values
        corr.covariance  = std::numeric_limits<double>::quiet_NaN();
        corr.correlation = std::numeric_limits<double>::quiet_NaN();

        // Cross product will be NaN or infinite if either x or y are.
        if( !finite(sum.xyw) || sum.count == 0 )
          continue;

        const double x_mean = sum.xw / sum.w;
        const double y_mean = sum.yw / sum.w;

        corr.covariance = sum.xyw / sum.w - x_mean * y_mean;

        if( fabs(corr.covariance) < 0.00000000001 )
          corr.covariance = 0.;

        if( sum.count < 3 )
          continue;

        const double x_variance  = sum.xxw / sum.w - x_mean * x_mean;
        const double y_variance  = sum.yyw / sum.w - y_mean * y_mean;
              double xy_variance = x_variance * y_variance;

        if( fabs(x_variance) < 0.00000000001 || fabs(y_variance) < 0.00000000001 )
          xy_variance = 0.;

        if( xy_variance > 0. )
          corr.correlation = corr.covariance / sqrt(xy_variance);

#if 0
  std::cout << x_variance << ", " << y_variance << ", " << xy_variance << "  ";
  std::cout << corr.covariance << ", " << corr.correlation << "  ";
  std::cout << std::endl;
#endif
     } // end of inner for loop.
   } // end of outer for loop.
} // end of update().

} // end of namespace SAGE
