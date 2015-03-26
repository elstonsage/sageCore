#ifndef HISTOGRAM_H
#define HISTOGRAM_H

// NOTE: This class needs to be completely re-written/structured.  The
//       current implementation is suitable for small-range integer samples. 
//       Large-range integers or real numbers will very likely result in
//       undesirable or erroneous behavoir.
 
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>
#include <cmath>
#include "sinfo.h"

namespace SAGE
{

//histograms group data points into bins and count the number the data
//points in each bin.
class Histogram
{
 public:

  //construct a histogram:
  // a.  max_bins is the maximum number of 'bins' or points in the histogram.
  // b.  reduce_factor is a factor which will convert the data points
  //     in the range 0 to max_bins-1
  // c.  Pre_condition  : max_bins > 0
  Histogram(size_t max_bins = 0, double reduce_factor = 1.0)
  {
    my_bound_flag = false;
    my_upper_bound = 0;
    my_lower_bound = 0;
    my_offset = 0;
    my_reduce_factor = reduce_factor;

    if(max_bins <= 0 )
    {
      my_data.resize(1);
      my_max_bins = 0;
      return;
    }

    my_max_bins = max_bins;
    my_data.resize(max_bins);
  };

  //construct a histogram given a range [min, max]
  // a.  max_bins is the maximum number of 'bins' or points in the histogram.
  // b.  Pre_condition  : max_bins > 0 and max > min
  Histogram(size_t max_bins, double min, double max, bool reduce = false)
  {
    my_bound_flag = true;
    my_reduce_factor = 1.0;
    my_offset = 0;

    if(max_bins <= 0 || min > max)     //treat it as no restrictions
    {
      my_max_bins = 0;
      my_bound_flag= false;
      my_lower_bound = my_upper_bound = 0;
      my_data.resize(1);
      return;
    }

    my_data.resize( max_bins );
    my_max_bins = max_bins;
    my_upper_bound = max;

    if(min < 0)
    {
      my_offset = -min;
      my_lower_bound = 0;
      my_upper_bound += my_offset;
    }
    else
    {
      //my_lower_bound = min;
      //my_offset = 0;

      // We still need offset to be set to -min when min >= 0
      // to add data into correct bin.  - bug fix by yjs 26-10-2001
      //
      my_offset = -min;
      my_lower_bound = 0;
      my_upper_bound += my_offset;
    }

    if(reduce)
      my_reduce_factor = (max-min)/max_bins;
  };

  //functions

  //Note: this operation can only be applied to the histograms with
  //      the same properties, that is:
  //       a.  same max_bins
  //       b.  same lower and upper bounds
  //       c.  same reduce factor
  Histogram&  operator +=( const Histogram& hg );
  Histogram   operator + ( const Histogram& hg ) const
  {
    Histogram temp(*this);
    return temp += hg;
  }

  double  binsize() const;
  bool    add(double x);

  //add 'points' to the bin 'index'. If the 'index' is invalid, return false
  bool    addPoints(size_t points, size_t index);

  //reset the histogram, max_bins and other conditions unaltered
  bool    clear();

  //return the count of items in bin with 'index'
  size_t     binCount(size_t index)  const ;

  // the maximum bin frequency and the bin that it occurs in
  std::pair<size_t,size_t> maxFrequency() const;

  Histogram&  operator += (size_t x)  { add(x); return *this;  }
  Histogram&  operator += (double x)  { add(x); return *this;  }

  size_t operator[] (int index)  const {  return my_data[index];  }

  bool    bounded()       const    { return my_bound_flag;          }
  size_t  max_bins()      const    { return my_max_bins;            }
  size_t  point_count()   const    { return my_datastat.count();    }
  size_t  total_point_count() const{ return point_count()
                                   + less_count() + more_count();     }
  size_t  less_count()    const    { return my_loutliers.count();     }
  size_t  more_count()    const    { return my_uoutliers.count();     }
  double  upper_bound()   const    { return my_upper_bound-my_offset; }
  double  lower_bound()   const    { return my_lower_bound-my_offset; }
  double  reduce_factor() const    { return my_reduce_factor;       }
  size_t  size()          const    { return my_data.size();         }
  double  mean()          const    { return my_datastat.mean();     }
  double  variance()      const    { return my_datastat.variance(); }

  std::pair<double,double> bin_range(size_t index) const;

  const SampleInfo&  stat_info()       const { return my_datastat;  }
  const SampleInfo   total_stat_info() const { return my_datastat
                                                    + my_loutliers
                                                    + my_uoutliers; }
  const SampleInfo&  lower_outliers()  const { return my_loutliers; }
  const SampleInfo&  upper_outliers()  const { return my_uoutliers; }

private:
  bool         my_bound_flag;
  size_t       my_max_bins;
  double       my_upper_bound;
  double       my_lower_bound;
  double       my_reduce_factor;
  double       my_offset; // = (0-my_lower_bound) if my_lower_bound<0
  std::vector<size_t> my_data;

  SampleInfo      my_datastat;
  SampleInfo      my_loutliers;   //samples that are outside of lower bound
  SampleInfo      my_uoutliers;   //samples that are outside of upper bound
};

//////////////////////////////
//                          //
//  inline implementations  //
//                          //
//////////////////////////////

inline double Histogram::binsize() const
{
 if(!bounded())   return 1.0/my_reduce_factor;
 else             return  ( upper_bound() - lower_bound() )/max_bins();
}

inline size_t  Histogram::binCount(size_t index)  const
{
 if(index >= my_data.size() )
  return std::numeric_limits<size_t>::quiet_NaN();
 else
  return my_data[index];
}

inline Histogram& Histogram::operator+=( const Histogram& hg )
{
  if(  max_bins() != hg.max_bins()|| reduce_factor() != hg.reduce_factor()
     ||lower_bound()!=hg.lower_bound() || upper_bound()!=hg.upper_bound() )
   return *this;

  my_datastat  = my_datastat  + hg.stat_info();
  my_loutliers = my_loutliers + hg.lower_outliers();
  my_uoutliers = my_uoutliers + hg.upper_outliers();

  //case: no max_bins, no boundaries
  bool flag = false;
  size_t  temp = size();

  if( max_bins() == 0 )
  {
    if( hg.size() > size() )
      flag=true;
    else
      temp = hg.size();
  }

  for(size_t i=0; i<temp; i++)
    addPoints(hg[i],i);

  if(flag)
    for(size_t j=temp; j<hg.size(); j++)
      addPoints(hg[j],j);

  return *this;
}

inline std::pair<size_t,size_t> Histogram::maxFrequency()  const
{
  size_t maxF=0;
  size_t index=(size_t)-1;

  for(size_t i=0; i<my_data.size(); i++)
    if(my_data[i] > maxF )
    {
      maxF=my_data[i];
      index=i;
    }
  return std::make_pair(index,maxF);
}

inline bool Histogram::add(double x)
{
  // no boundaries, no max_bins
  x += my_offset;

  if(max_bins() ==  0)
  {
    if(floor(x) >= my_data.size()) my_data.resize( (size_t)floor(x)+1);
    my_data[ (size_t)floor(x) ] += 1;
    my_datastat.add(x);
    return true;
  }

  //bounded
  if( bounded() )
  {
    if(x > my_upper_bound)  { my_uoutliers.add(x); return false;}
    if(x < my_lower_bound)  { my_loutliers.add(x); return false;}

    size_t bin_index = (size_t)floor(x/binsize());

    my_data[ bin_index ] += 1;
    my_datastat.add(x);
    return true;
  }

  //not bounded, but with max_bins and reduce_factor
  if( floor(x*my_reduce_factor) >= max_bins()) return false;
  my_data[ (size_t)floor(x*my_reduce_factor) ] += 1;
  my_datastat.add(x);

  return true;
}

inline std::pair<double,double> Histogram::bin_range(size_t index) const
{
  double qNaN = std::numeric_limits<double>::quiet_NaN();
  double eps  = std::numeric_limits<double>::epsilon();

  if(index >= size())
    return std::make_pair(qNaN,qNaN);

  // no boundaries, no max_bins
  if(max_bins() ==  0)
    return std::make_pair(index, index+1.0-eps);

  //bounded
  if( bounded() )
  {
    double r   = (my_upper_bound-my_lower_bound)/(double)max_bins();
    return std::make_pair( index*r, (index+1)*r-eps );
  }

  //not bounded, but with max_bins and reduce_factor
  return std::make_pair( index/my_reduce_factor, (index+1)/my_reduce_factor-eps);
}

inline  bool Histogram::addPoints(size_t points, size_t index)
{
  if( points==0)  return true;

  if(bounded() && ( index < lower_bound() || index > upper_bound() ) )
    return false;

  if(max_bins() &&  index >= max_bins() )
    return false;

  if(! add(index) )  return false;

  for(size_t i=1; i<points; i++)
    add(index);

  return true;
}

inline bool Histogram::clear()
{
  my_data.clear();
  my_datastat.clear();
  my_datastat.clear();
  my_loutliers.clear();
  my_uoutliers.clear();
  return true;
}

}

#endif
