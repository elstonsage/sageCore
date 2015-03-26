#include <cmath>

//lint --e{1727}
//lint --e{1529}

static const double MAX_DOUBLE_SCALE = sqrt(std::numeric_limits<double>::max());
static const double MIN_DOUBLE_SCALE = sqrt(std::numeric_limits<double>::min());

static const double LOG_MAX_DOUBLE_SCALE = log(sqrt(std::numeric_limits<double>::max()));
static const double LOG_MIN_DOUBLE_SCALE = log(sqrt(std::numeric_limits<double>::min()));

static const double MAX_ADDITIVE_LOGS    = log(10.0 / sqrt(std::numeric_limits<double>::epsilon()));

inline
log_double::log_double()
  : is_on_double_scale(false),
    doub(-std::numeric_limits<double>::infinity())
{ }

inline
log_double::log_double(const log_double& d)
  : is_on_double_scale(d.is_on_double_scale),
    doub(d.doub)
{ }

inline
log_double::log_double(double d)
  : is_on_double_scale(true),
    doub(d)
{ }

inline
log_double& log_double::operator=(const log_double& d)
{
  is_on_double_scale = d.is_on_double_scale;

  doub = d.doub;

  return *this;
}

inline
log_double& log_double::operator=(double d)
{
  is_on_double_scale = true;

  doub = d;

  convert();

  return *this;
}

inline
double log_double::get_double() const
{
  if(is_on_double_scale) return doub;
  else                   return exp(doub);
}

inline
double log_double::get_log() const
{
  if(is_on_double_scale) return log(doub);
  else                   return doub;
}

inline
log_double log_double::operator*(const log_double& d) const
{
  log_double temp(*this);

  temp *= d;

  return temp;
}

inline
log_double log_double::operator*(double d) const
{
  log_double temp(*this);

  temp *= d;

  return temp;
}

inline
log_double& log_double::operator*=(const log_double& d)
{
  if(is_on_double_scale)
  {
    if(d.is_on_double_scale)
    {
      doub = doub * d.doub;
    }
    else
    {
      is_on_double_scale = false;

      doub = log(doub) + d.doub;
    }
  }
  else
  {
    if(d.is_on_double_scale)
    {
      doub += log(d.doub);
    }
    else
    {
      doub += d.doub;
    }
  }

  convert();

  return *this;
}

inline
log_double& log_double::operator*=(double d)
{
  *this *= log_double(d);

  return *this;
}

inline
log_double log_double::operator/(const log_double& d) const
{
  log_double temp(*this);

  temp /= d;

  return temp;
}

inline
log_double log_double::operator/(double d) const
{
  log_double temp(*this);

  temp /= d;

  return temp;
}

inline
log_double& log_double::operator/=(const log_double& d)
{
  if(is_on_double_scale)
  {
    if(d.is_on_double_scale)
    {
      doub /= d.doub;
    }
    else
    {
      is_on_double_scale = false;

      doub = log(doub) - d.doub;
    }
  }
  else
  {
    if(d.is_on_double_scale)
    {
      doub -= log(d.doub);
    }
    else
    {
      doub -= d.doub;
    }
  }

  convert();

  return *this;
}

inline
log_double& log_double::operator/=(double d)
{
  *this /= log_double(d);

  return *this;
}

inline
log_double log_double::operator+(const log_double& d) const
{
  log_double temp(*this);

  temp += d;

  return temp;
}

inline
log_double log_double::operator+(double d) const
{
  log_double temp(*this);

  temp += log_double(d);

  return temp;
}

inline
log_double& log_double::operator+=(const log_double& d)
{
  if(is_on_double_scale)
  {
    if(d.is_on_double_scale)
    {
      doub += d.doub;
    }
    else
    {
      is_on_double_scale = false;
      doub = add_logs(log(doub),d.doub);
    }
  }
  else
  {
    if(d.is_on_double_scale)
    {
      doub = add_logs(doub, log(d.doub));
    }
    else
    {
      doub = add_logs(doub,d.doub);
    }
  }

  convert();

  return *this;
}

inline
log_double& log_double::operator+=(double d)
{
  *this += log_double(d);

  return *this;
}


inline
log_double log_double::operator-(const log_double& d) const
{
  log_double temp(*this);

  temp -= d;

  return temp;
}

inline
log_double log_double::operator-(double d) const
{
  log_double temp(*this);

  temp -= log_double(d);

  return temp;
}

inline
log_double& log_double::operator-=(const log_double& d)
{
  if(is_on_double_scale)
  {
    if(d.is_on_double_scale)
    {
      doub -= d.doub;
    }
    else
    {
      is_on_double_scale = false;
      doub = subtract_logs(log(doub),d.doub);
    }
  }
  else
  {
    if(d.is_on_double_scale)
    {
      doub = subtract_logs(doub, log(doub));
    }
    else
    {
      doub = subtract_logs(doub,d.doub);
    }
  }

  convert();

  return *this;
  
}

inline
log_double& log_double::operator-=(double d)
{
  *this -= log_double(d);

  return *this;
}


inline
log_double& log_double::pow(double d)
{
  if(is_on_double_scale)
  {
    is_on_double_scale = false;
    doub = log(doub) * d;
  }
  else
    doub *= d;

  convert();

  return *this;
}

inline
bool log_double::compare_equals(const log_double& rhs, const log_double& epsilon) const
{
  return *this >= (rhs - epsilon) && *this <= (rhs + epsilon);
}

inline
bool log_double::compare_equals(double rhs, double epsilon) const
{
  return compare_equals(log_double(rhs), log_double(epsilon));
}

inline
bool log_double::compare_not_equals(const log_double& rhs, const log_double& epsilon) const
{
  return !compare_equals(rhs, epsilon);
}

inline
bool log_double::compare_not_equals(double rhs, double epsilon) const
{
  return !compare_equals(rhs, epsilon);
}

inline
bool log_double::operator==(const log_double& rhs) const
{
  //lint --e{777}
  return is_on_double_scale == rhs.is_on_double_scale && doub == rhs.doub;
}

inline
bool log_double::operator!=(const log_double& rhs) const
{
  //lint --e{777}
  return is_on_double_scale != rhs.is_on_double_scale || doub != rhs.doub;
}

inline
bool log_double::operator< (const log_double& d) const
{
  return this->get_log() < d.get_log();
}

inline
bool log_double::operator< (double d) const
{
  return *this < log_double(d);
}

inline
bool log_double::operator<=(const log_double& d) const
{
  return this->get_log() <= d.get_log();
}

inline
bool log_double::operator<=(double d) const
{
  return *this <= log_double(d);
}

inline
bool log_double::operator> (const log_double& d) const
{
  return this->get_log() > d.get_log();
}

inline
bool log_double::operator> (double d) const
{
  return *this > log_double(d);
}

inline
bool log_double::operator>=(const log_double& d) const
{
  return this->get_log() >= d.get_log();
}

inline
bool log_double::operator>=(double d) const
{
  return *this >= log_double(d);
}

inline 
void log_double::convert()
{
  if(is_on_double_scale)
  {
    if(doub < MIN_DOUBLE_SCALE || doub > MAX_DOUBLE_SCALE)
    {
      doub = log(doub);
      is_on_double_scale = false;
    }
  }
  else
  {
    if(doub > LOG_MIN_DOUBLE_SCALE && doub < LOG_MAX_DOUBLE_SCALE)
    {
      doub = exp(doub);
      is_on_double_scale = true;
    }
  }
}

inline
double log_double::add_logs(double d1, double d2) const
{
  double difference = d1 - d2;

  if(difference < 0.0) difference = -difference;

  // If the difference is not finite, that means that one of the values is
  // not finite (NAN or infinity).
  if(!finite(difference))
  {
    if(SAGE::isnan(d1) || SAGE::isnan(d2))
    {
      d1 = std::numeric_limits<double>::quiet_NaN();
    }
    else if(d2 > d1)
    {
      d1 = d2;
    }

    return d1;
  }

  // If the difference is so large it wouldn't show with the double, we just
  // return the bigger value.
  if(difference / 2.0 > MAX_ADDITIVE_LOGS)
  {
    if(d2 > d1) return d2;
    else        return d1;
  }

  // Otherwise, we can actually add the values, but we do this by
  // centering at 0.  This means that the values will be +/- difference / 2,
  // and the values on the double scale x and 1/x.  We thus only have to exp
  // once.

  double standard = (d1 + d2) / 2;

  double diff_exp = exp(difference / 2);

  d1 = log(diff_exp + 1.0 / diff_exp) + standard;

  return d1;
}

inline
double log_double::subtract_logs(double d1, double d2) const
{
  if(d2 > d1)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double difference = d1 - d2;

  // If the difference is NaN, that means that one of the values is
  // not finite.
  if(!finite(difference))
  {
    if(SAGE::isnan(d1) || SAGE::isnan(d2))
    {
      d1 = std::numeric_limits<double>::quiet_NaN();
    }

    return d1;
  }

  // If the difference is so large it wouldn't show with the double, we
  // return the first value.
  if(difference/2.0 > MAX_ADDITIVE_LOGS)
  {
    return d1;
  }

  // Otherwise, we can actually add the values, but we do this by
  // centering at 0.  This means that the values will be +/- difference / 2,
  // and the values on the double scale x and 1/x.  We thus only have to exp
  // once.

  double standard = (d1 + d2) / 2;

  double diff_exp = exp(difference / 2);

  d1 = log(diff_exp - 1.0 / diff_exp) + standard;

  return d1;
}
