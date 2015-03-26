#include <iostream>
#include <iomanip>
#include <math.h>
#include <limits>
#include "numerics/log_double.h"
#include "numerics/constants.h"

namespace SAGE {
namespace NUMERICS
{

/// The normal_pdf function calculates the normal point density function
/// given a mean and a variance.  The equation for this is:
///
/// \f[ \frac{1} {\sqrt{ 2 \pi w } } exp \left(  - \frac{ z^2 } { 2 w } \right) \f]
///
/// The assumption is that mean and var are finite real numbers, and
/// furthermore that var is > 0.  If these criteria are not met, then the
/// behavior is undefined.
///
/// \param mean The mean for which the normal PD is desired.
/// \param var  The variance for which the normal PD is desired.
inline 
double
normal_pdf(double mean, double var)
{
    // Calculate the sqrt component of variance.

    double sqrt_component = 1.0 / sqrt(2.0 * numeric_constants<double>::pi() * var);

    // Calculate the exponential component of mean and variance.

    double exp_component = exp( -0.5 * mean * mean / var );

    // The PDF is the product of the two components

    return sqrt_component * exp_component;
}

/// The log_normal_pdf function implments the normal_pdf() function,
/// but returns the calculated value as a log_double.
///
/// This is useful if you're running into precision problems with normal_pdf.
///
/// \param mean The mean for which the normal PD is desired.
/// \param var  The variance for which the normal PD is desired.
inline 
log_double
log_normal_pdf(double mean, double var)
{
    // Calculate the sqrt component of variance.

    double sqrt_component = 1.0 / sqrt(2.0 * numeric_constants<double>::pi() * var);

    // Calculate the exponential component of mean and variance.

    double exp_component = exp( -0.5 * mean * mean / var );

    // The PDF is the product of the two components
    return log_double(sqrt_component * exp_component);
}

/// The normal_pdf function calculates the normal point density function
/// given a mean and a variance.
///
/// The assumption is that mean and var are finite real numbers, and
/// furthermore that var is > 0.  If these criteria are not met, then the
/// behavior is undefined.
///
/// \param z The z score
/// \param mean The mean
/// \param var The variance
inline double
normal_pdf(double z, double mean, double var)
{
    // Calculate the sqrt component of variance.

    double sqrt_component = 1.0 / sqrt(2.0 * numeric_constants<double>::pi() * var);

    // Calculate the exponential component of mean and variance.

    double exp_component = exp(-0.5 * ((z - mean) * (z - mean)) / var);

    // The PDF is the product of the two components

    return sqrt_component * exp_component;
}

} // End NUMERICS namespace
} // End SAGE namespace
