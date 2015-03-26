#ifndef ASSOC_RESIDUALS_H
#define ASSOC_RESIDUALS_H

//============================================================================
// File:      Residuals.h
//
// Author:    Dan Baechle
//
// History:   6/6/11 - created.                                   djb
//
// Notes:     Represents a set of residuals in the Assoc regression
//            model.
//
//
// Copyright (c) 2011 R.C. Elston
// All Rights Reserved
//============================================================================


#include <vector>
#include <cmath>
#include "numerics/sinfo.h"

using std::vector;
using std::isnan;

namespace SAGE  
{

namespace ASSOC 
{

class Residuals
{
  public:
    Residuals();
    Residuals(const vector<double>& raw_residuals, double rescaling_factor);

    size_t  size() const;
    
    double  operator[](size_t index) const;
    double  standardized_residual(size_t index) const;
    
    double  mean() const;
    double  std_dev() const;
    double  variance() const;
    double  skewness() const;
    double  kurtosis() const;
    
  private:
    vector<double>  my_data;
    SampleInfo  my_stats;
};

#include "assoc/Residuals.ipp"

}
}
  
#endif


