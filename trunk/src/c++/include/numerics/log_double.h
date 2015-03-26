#ifndef LOG_DOUBLE_H
#define LOG_DOUBLE_H

#include "numerics/isnan.h"
#include "globals/config.h"

#include <math.h>
#include <limits>
#include <iostream>

/** log_double stores a double on the log scale
 *  This class allows storing of doubles on a logarithmic scale.  It is assumed
 *  that the numbers are regarded as normal doubles, in that the operations * +, etc.
 *  are applied to the underlying doubles.  Under this, * becomes addition of the logs,
 *  and + becomes complicated.
 */
class log_double
{
  friend std::ostream& operator<<(std::ostream& o, const log_double& d);


  //lint --e{1739}

  public:

             log_double();
             log_double(const log_double&);
    explicit log_double(double);

    log_double& operator=(const log_double&);
    log_double& operator=(double);
    
    double get_double() const;
    double get_log   () const;

    // Operators
    
    log_double operator*(const log_double&) const;
    log_double operator*(double)            const;

    log_double& operator*=(const log_double&);
    log_double& operator*=(double);
    
    log_double operator/(const log_double&) const;
    log_double operator/(double)            const;

    log_double& operator/=(const log_double&);
    log_double& operator/=(double);
  
    log_double operator+(const log_double&) const;
    log_double operator+(double)            const;

    log_double& operator+=(const log_double&);
    log_double& operator+=(double);
    
    log_double operator-(const log_double&) const;
    log_double operator-(double)            const;

    log_double& operator-=(const log_double&);
    log_double& operator-=(double);

    // Power function:
    
    log_double& pow(double);

    // Comparisons
    
    bool compare_equals(const log_double& rhs, const log_double& epsilon) const;
    bool compare_equals(double            rhs,            double epsilon) const;

    bool compare_not_equals(const log_double& rhs, const log_double& epsilon) const;
    bool compare_not_equals(double            rhs,            double epsilon) const;

    bool operator==(const log_double&) const;
    bool operator!=(const log_double&) const;

    bool operator< (const log_double&) const;
    bool operator< (double)            const;

    bool operator<=(const log_double&) const;
    bool operator<=(double)            const;

    bool operator> (const log_double&) const;
    bool operator> (double)            const;

    bool operator>=(const log_double&) const;
    bool operator>=(double)            const;

  private:

    void   convert      ();
    double add_logs     (double d1, double d2) const;
    double subtract_logs(double d1, double d2) const;

    bool   is_on_double_scale;

    double doub;
};

inline log_double operator*(double d, const log_double& d2)
{
  return d2 * d;
}

inline std::ostream& operator<<(std::ostream& o, const log_double& d)
{
  o << d.get_double() << "//" << d.get_log();

  return o;
}

#include "numerics/log_double.ipp"

#endif
