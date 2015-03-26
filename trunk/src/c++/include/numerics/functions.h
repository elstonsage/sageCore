#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// include the CEPHES numerics
#include <stdlib.h>
#include "numerics/cephes.h"

namespace SAGE
{


///
/// Given a double v, returns the number of decimal places to the left of
/// the decimal point (excluding the sign).
inline int getLeftLength(double v) { int int_v = abs((int)v); return int_v > 0 ? (int) log10((double)int_v) + 1: 1; }

inline double binomial_cdf(int k, int n, double x)     { return bdtr(k,n,x); }
inline double neg_binomial_cdf(int k, int n, double x) { return nbdtr(k,n,x); }

// Choice is good.

inline double chi_square_cdf(int df, double x)     { return chdtr(df,x);  }
inline double inv_chi_square_cdf(int df, double x) { return chdtri(df,x); }

inline double F_cdf(int df1, int df2, double x) { return fdtr(df1,df2,x); }
inline double inv_F_cdf(int df1, int df2, double p) { return fdtri(df1,df2,p);}

inline double log_gamma(double x)                      { return lgam(x); }
inline double incomplete_gamma(double a, double x)     { return igam(a,x); }
inline double inv_incomplete_gamma(double a, double p) { return igami(a,p); }
inline double gamma_cdf(double a, double b, double x)  { return gdtr(a,b,x); }

// Default assumption is mean=0, sdev=1
inline double normal_cdf(double x)     { return ndtr(x); }
inline double inv_normal_cdf(double p) { return ndtri(p); }

// Convenience routines for other values of the mean and sdev
inline double normal_cdf(double m, double s, double x) 
{
  return s > 0 ? ndtr((x-m)/s) : (x<m ? 0 : (x>m ? 1 : 0.5));
}

inline double inv_normal_cdf(double m, double s, double p) 
{ return m + s*ndtri(p); }

inline double poisson_cdf(int k, double x) { return pdtr(k,x); }
inline double inv_poisson_cdf(int k, double p) { return pdtri(k,p); }

// You can use it either way...
inline double students_cdf(int df, double t)     { return stdtr(df,t); }
inline double T_cdf(int df, double t)            { return stdtr(df,t); }
inline double inv_students_cdf(int df, double p) { return stdtri(df,p); }
inline double inv_T_cdf(int df, double p)        { return stdtri(df,p); }

/* probability density from -t to t */
inline double T_mean_test(int k, double t)       { return stdtr_mean_test(k, t); }

/* t for range -t to t for probability density p */ 
inline double T_mean_test_inv(int k, double p)   { return stdtr_mean_test_inv(k, p); }
}

#endif 
