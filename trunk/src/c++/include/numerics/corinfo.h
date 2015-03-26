#ifndef CORRELATIONINFO_H
#define CORRELATIONINFO_H

//****************************************************************************
//* File:      corinfo.h                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song & Kevin Jacobs      			     *
//*					                                     *
//* History:   Version 1.0 						     *
//*									     *
//* Notes:     This header file defines correlation information class to be  *
//*            used to calculate the correlations of (xi,xj) pairs of        *
//*            multivariates (x1, x2, x3, ..., xN).                          *
//*                       						     *
//* Copyright (c) 1999 R.C. Elston 					     *
//*   All Rights Reserved 						     *
//****************************************************************************

#include <limits>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <vector>
#include "numerics/kahan.h"
#include "numerics/sinfo.h"
#include "numerics/trimatrix.h"

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SimpleCorrelationInfo                                        ~
// ~                                                                         ~
// ~ Purpose:   Compute the following statistics for a given sample          ~
// ~              (x1, x2)                                                   ~
// ~             * covariance  - covariance of (x1,x2)                       ~
// ~             * correlation - correlation of (xi, xj)                     ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SimpleCorrelationInfo
{
  public:

    typedef KahanAdder<double>                                sum_type;

    SimpleCorrelationInfo();

    double            mean1()                          const;
    double            variance1()                      const;
    const SampleInfo& info1()                          const;

    double            mean2()                          const;
    double            variance2()                      const;
    const SampleInfo& info2()                          const;

    double            correlation()                    const;
    double            covariance()                     const;

    size_t            count()                          const;

    void              add(double i, double j);
    
    void              clear();

    SimpleCorrelationInfo&  operator+=(const SimpleCorrelationInfo &c);
    SimpleCorrelationInfo   operator+ (const SimpleCorrelationInfo &c)               const;

  private:


    struct correlation_computation_type
    {
      correlation_computation_type();
      correlation_computation_type& operator+=(const correlation_computation_type &c);
      correlation_computation_type  operator+ (const correlation_computation_type &c) const;

      void add(double a1, double a2, double wi = 1.0);

      sum_type   x;
      sum_type   y;
      sum_type   xx;
      sum_type   yy;
      sum_type   xy;
    };

    struct correlation_type
    {
      correlation_type();

      double     covariance;
      double     correlation;
    };

    size_t                                        my_count;

    SampleInfo                                    my_info1;
    SampleInfo                                    my_info2;

    correlation_computation_type                  my_sums;

    // Mutable cache of computed values and when they were last computed.
    mutable correlation_type                      my_correlation;
    mutable size_t                                my_last_count;

    // Recompute correlations and covariances
    void    update()                              const;
    void    do_update()                           const;

}; // end of class definition

// ---------------------------------------------------------------------------
// Implementation of SimpleCorrelationInfo::correlation_computation_type
// ---------------------------------------------------------------------------

inline
SimpleCorrelationInfo::correlation_computation_type::
correlation_computation_type() : x(0), y(0), xx(0), yy(0), xy(0)
{ }

inline
SimpleCorrelationInfo::correlation_computation_type&
SimpleCorrelationInfo::correlation_computation_type::
operator+=(const correlation_computation_type &c)
{
    x      +=  c.x;
    y      +=  c.y;
    xx     +=  c.xx;
    yy     +=  c.yy;
    xy     +=  c.xy;

    return *this;
}

inline
SimpleCorrelationInfo::correlation_computation_type
SimpleCorrelationInfo::correlation_computation_type::
operator+ (const correlation_computation_type &c) const
{
    correlation_computation_type temp(*this);
    return temp += c;
}

inline
void
SimpleCorrelationInfo::correlation_computation_type::
add(double a1, double a2, double wi)
{
    if( !SAGE::isnan(a1) && !SAGE::isnan(a2) )
    {
      x       +=  a1;
      y       +=  a2;
      xx      +=  a1 * a1;
      yy      +=  a2 * a2;
      xy      +=  a1 * a2;
    }
}

inline
SimpleCorrelationInfo::correlation_type::correlation_type() 
              : covariance(std::numeric_limits<double>::quiet_NaN()),
                correlation(std::numeric_limits<double>::quiet_NaN())
{ }

// ---------------------------------------------------------------------------
// Implementation of SimpleCorrelationInfo
// ---------------------------------------------------------------------------

inline
SimpleCorrelationInfo::SimpleCorrelationInfo()
{
    my_count = 0;
    my_last_count = 0;
}

inline
void
SimpleCorrelationInfo::clear()
{
    my_info1.clear();
    my_info2.clear();

    my_count = 0;
    my_sums = correlation_computation_type();
    my_correlation = correlation_type();
    my_last_count = 0;
}

inline
const SampleInfo&
SimpleCorrelationInfo::info1() const
{
    return my_info1;
}

inline
const SampleInfo&
SimpleCorrelationInfo::info2() const
{
    return my_info2;
}

inline
double
SimpleCorrelationInfo::mean1() const
{
    return my_info1.mean();
}

inline
double
SimpleCorrelationInfo::mean2() const
{
    return my_info2.mean();
}

inline
double
SimpleCorrelationInfo::variance1() const
{
  return my_info1.variance();
}

inline
double
SimpleCorrelationInfo::variance2() const
{
  return my_info2.variance();
}

inline
double
SimpleCorrelationInfo::correlation() const
{
    update();
    return my_correlation.correlation;
}

inline
double
SimpleCorrelationInfo::covariance() const
{
    update();

    return my_correlation.covariance;
}

inline
size_t
SimpleCorrelationInfo::count() const
{
    return my_count;
}

inline
SimpleCorrelationInfo
SimpleCorrelationInfo::operator+(const SimpleCorrelationInfo &c) const
{
    SimpleCorrelationInfo temp(*this);
    return temp += c;
}

inline void
SimpleCorrelationInfo::update() const
{
    if(my_count == my_last_count)
      return;

    do_update();
}

// end of SimpleCorrelationInfo Implementation

class CorrelationInfo
{
  public:

    typedef KahanAdder<double>                                sum_type;

    CorrelationInfo();
    CorrelationInfo(size_t variates);

    double            mean(size_t i)                          const;
    double            variance(size_t i)                      const;

    double            correlation(size_t i)                   const;
    double            covariance(size_t i)                    const;

    double            correlation(size_t i, size_t j)         const;
    double            covariance(size_t i, size_t j)  	      const;

    sum_type          sum_weight(size_t i)                    const;
    sum_type          sum_weight_square(size_t i)             const;

    sum_type          sum_weight(size_t i, size_t j)          const;
    sum_type          sum_weight_square(size_t i, size_t j)   const;

    size_t            size()                                  const;
    size_t            count()                                 const;
    size_t            count(size_t i)                         const;
    size_t            count(size_t i, size_t j)               const;

    void              add(const std::vector<double>& data, double weight = 1.0);
    void              add(size_t i, double di, size_t j, double dj, double weight = 1.0);
    void              add(const std::vector<double>& data, const TriangleMatrix<double>& weight);
    
    void              clear();
    void              resize(size_t variates);

    CorrelationInfo&  operator+=(const CorrelationInfo &c);
    CorrelationInfo   operator+ (const CorrelationInfo &c)               const;

    // Can be used only for unweighted corinfo.
    CorrelationInfo&  operator+=(const std::vector<double>& data);
    CorrelationInfo   operator+ (const std::vector<double>& data)        const;

    void    correlation(TriangleMatrix< KahanAdder<double> >& t)         const;
    void    covariance (TriangleMatrix< KahanAdder<double> >& t)         const;
    
    void    correlation(TriangleMatrix< double >& t)                     const;
    void    covariance (TriangleMatrix< double >& t)                     const;
    
    // [kbj] I don't like this!
    //       -- we should expose correlation_computation_type if we
    //       -- want this kind of functionality!
    
    void    sum_weight       (TriangleMatrix< KahanAdder<double> >& t)   const;
    void    sum_weight_square(TriangleMatrix< KahanAdder<double> >& t)   const;

    void    sum_weight       (TriangleMatrix<double>& t)                 const;
    void    sum_weight_square(TriangleMatrix<double>& t)                 const;
    
  private:


    struct correlation_computation_type
    {
      correlation_computation_type();

      correlation_computation_type& operator+=(const correlation_computation_type &c);

      correlation_computation_type  operator+ (const correlation_computation_type &c) const;

      void add(double a1, double a2, double wi = 1.0);

      size_t     count;
      
      sum_type   xw;
      sum_type   yw;
      sum_type   xxw;
      sum_type   yyw;
      sum_type   xyw;
      sum_type   w;
      sum_type   ww;
    };

    struct correlation_type
    {
      correlation_type();

      double     covariance;
      double     correlation;
    };

    size_t                                        my_count;

    std::vector<SampleInfo>                       my_variates;
    TriangleMatrix<correlation_computation_type>  my_sums;

    // Mutable cache of computed values and when they were last computed.
    mutable TriangleMatrix<correlation_type>      my_correlations;
    mutable size_t                                my_last_count;

    // Recompute correlations and covariances
    void    update()                              const;
    void    do_update()                           const;

}; // end of class definition

// ---------------------------------------------------------------------------
// Implementation of CorrelationInfo::correlation_computation_type
// ---------------------------------------------------------------------------

inline
CorrelationInfo::correlation_computation_type::
correlation_computation_type() : count(0), xw(0), xxw(0), yyw(0), xyw(0), w(0), ww(0)
{ }

inline
CorrelationInfo::correlation_computation_type&
CorrelationInfo::correlation_computation_type::
operator+=(const correlation_computation_type& c)
{
    count   +=  c.count;
    xw      +=  c.xw;
    yw      +=  c.yw;
    xxw     +=  c.xxw;
    yyw     +=  c.yyw;
    xyw     +=  c.xyw;
    w       +=  c.w;
    ww      +=  c.ww;

    return *this;
}

inline
CorrelationInfo::correlation_computation_type
CorrelationInfo::correlation_computation_type::
operator+ (const correlation_computation_type &c) const
{
    correlation_computation_type temp(*this);
    return temp += c;
}

inline
void
CorrelationInfo::correlation_computation_type::
add(double a1, double a2, double wi)
{
    if( !SAGE::isnan(a1) && !SAGE::isnan(a2) )
    {
      count   +=  1;
      xw      +=  a1 * wi;
      yw      +=  a2 * wi;
      xxw     +=  a1 * a1 * wi;
      yyw     +=  a2 * a2 * wi;
      xyw     +=  a1 * a2 * wi;
      w       +=  wi;
      ww      +=  wi * wi;
    }
}

inline
CorrelationInfo::correlation_type::correlation_type() 
              : covariance(std::numeric_limits<double>::quiet_NaN()),
                correlation(std::numeric_limits<double>::quiet_NaN())
{ }

// ---------------------------------------------------------------------------
// Implementation of CorrelationInfo
// ---------------------------------------------------------------------------

inline
CorrelationInfo::CorrelationInfo()
{
    my_count = 0;
    my_last_count = 0;
}

inline
CorrelationInfo::CorrelationInfo(size_t number_of_variates)
{
    resize(number_of_variates);
}

inline
void
CorrelationInfo::clear()
{
    resize(0);
    resize(my_variates.size());
}

inline
double
CorrelationInfo::mean(size_t i) const
{
    return my_variates[i].mean();
}

inline
double
CorrelationInfo::variance(size_t i) const
{
  return covariance(i,i);
}

inline
double
CorrelationInfo::correlation(size_t i) const
{
    if( i >= my_correlations.linear_size() )
      return std::numeric_limits<double>::quiet_NaN();

    update();

    return my_correlations[i].correlation;
}

inline
double
CorrelationInfo::covariance(size_t i) const
{
    if( i >= my_correlations.linear_size() )
      return std::numeric_limits<double>::quiet_NaN();

    update();

    return my_correlations[i].covariance;
}

inline
double
CorrelationInfo::correlation(size_t i, size_t j) const
{
    if( std::max(i,j) >= size() )
      return std::numeric_limits<double>::quiet_NaN();

    if( i == j )
      return 1.0;

    update();

    return my_correlations(i, j).correlation;
}

inline
double
CorrelationInfo::covariance(size_t i, size_t j) const
{
    if( std::max(i,j) >= size() )
      return std::numeric_limits<double>::quiet_NaN();

    if( i == j )
      return my_variates[i].variance();

    update();

    return my_correlations(i, j).covariance;
}

inline
CorrelationInfo::sum_type
CorrelationInfo::sum_weight(size_t i) const
{
    if( i >= my_sums.linear_size() )
      return sum_type(0);

    return my_sums[i].w;
}

inline
CorrelationInfo::sum_type
CorrelationInfo::sum_weight_square(size_t i) const
{
    if( i >= my_sums.linear_size() )
      return sum_type(0);

    return my_sums[i].ww;
}

inline
CorrelationInfo::sum_type
CorrelationInfo::sum_weight(size_t i, size_t j) const
{
    if( std::max(i,j) >= size() )
      return sum_type(0);

    return my_sums(i, j).w;
}

inline
CorrelationInfo::sum_type
CorrelationInfo::sum_weight_square(size_t i, size_t j) const
{
    if( std::max(i,j) >= size() )
      return sum_type(0);

    return my_sums(i, j).ww;
}

inline
size_t
CorrelationInfo::size() const
{
    return my_variates.size();
}

inline
size_t
CorrelationInfo::count() const
{
    return my_count;
}

inline
size_t
CorrelationInfo::count(size_t i) const
{
    if( i >= size() )
      return 0;

    return my_variates[i].count();
}

inline
size_t
CorrelationInfo::count(size_t i, size_t j) const
{
    if( i == j )
      return count(i);

    if( i >= size() )
      return 0;

    return my_sums(i, j).count;
}


inline
CorrelationInfo
CorrelationInfo::operator+(const CorrelationInfo &c) const
{
    CorrelationInfo temp(*this);
    return temp += c;
}

inline
CorrelationInfo&
CorrelationInfo::operator+=(const std::vector<double>& data)
{
  add(data);
  return *this;
}

inline
CorrelationInfo
CorrelationInfo::operator+(const std::vector<double>& data) const
{
    CorrelationInfo temp(*this);
    return temp += data;
}

inline
void
CorrelationInfo::correlation(TriangleMatrix< KahanAdder<double> >& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    update();

    for( size_t i = 0; i < my_correlations.linear_size(); ++i )
        t[i] = my_correlations[i].correlation;
}           

inline
void
CorrelationInfo::covariance(TriangleMatrix< KahanAdder<double> >& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    update();

    for( size_t i = 0; i < my_correlations.linear_size(); ++i )
        t[i] = my_correlations[i].covariance;
}           

inline
void
CorrelationInfo::correlation(TriangleMatrix< double >& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    update();

    for( size_t i = 0; i < my_correlations.linear_size(); ++i )
        t[i] = my_correlations[i].correlation;
}           

inline
void
CorrelationInfo::covariance(TriangleMatrix< double >& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    update();

    for( size_t i = 0; i < my_correlations.linear_size(); ++i )
        t[i] = my_correlations[i].covariance;
}           

inline
void
CorrelationInfo::sum_weight(TriangleMatrix<double>& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    for( size_t i = 0; i < my_sums.linear_size(); ++i )
        t[i] = my_sums[i].w;
}           

inline
void
CorrelationInfo::sum_weight_square(TriangleMatrix<double>& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    for( size_t i = 0; i < my_sums.linear_size(); ++i )
        t[i] = my_sums[i].ww;
}           

inline
void
CorrelationInfo::sum_weight(TriangleMatrix< KahanAdder<double> >& t) const
{
    if( t.size() != size() )
      t.resize(size());

    for( size_t i = 0; i < my_sums.linear_size(); ++i )
        t[i] = my_sums[i].w;
}           

inline
void
CorrelationInfo::sum_weight_square(TriangleMatrix< KahanAdder<double> >& t) const
{
    if( t.size() != size() )
      t.resize(size());
      
    for( size_t i = 0; i < my_sums.linear_size(); ++i )
        t[i] = my_sums[i].ww;
}           

inline void
CorrelationInfo::update() const
{
    if(my_count == my_last_count)
      return;

    do_update();
}

// end of CorrelationInfo Implementation

} // end of namespace SAGE

#endif
