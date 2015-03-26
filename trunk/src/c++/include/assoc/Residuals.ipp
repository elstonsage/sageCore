//============================================================================
// File:      Residuals.ipp
//
// Author:    Dan Baechle
//
// History:   6/6/11 - created.                                   djb
//
// Notes:    
//
//
// Copyright (c) 2011 R.C. Elston
// All Rights Reserved
//============================================================================


//============================================================================
// IMPLEMENTATION:  Residuals       
//============================================================================
//
inline size_t
Residuals::size() const
{
  return  my_data.size();
}
        

inline  double
Residuals::operator[](size_t index) const
{
  return  my_data[index];
}


inline  double
Residuals::standardized_residual(size_t index) const
{
  return  (my_data[index] - mean()) / std_dev();
}
                
inline  double  
Residuals::mean() const
{
  return  my_stats.mean();
}

inline double  
Residuals::std_dev() const
{
  return  my_stats.standard_deviation();
}
    
inline double  
Residuals::variance() const
{
  return  my_stats.variance();
}

inline double  
Residuals::skewness() const
{
  return  my_stats.skewness();
} 
 
inline double  
Residuals::kurtosis() const
{
  return  my_stats.kurtosis();
}

