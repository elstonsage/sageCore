//============================================================================
// File:      Residuals.cpp
//
// Author:    Dan Baechle
//
// History:   created 6/6/11                               djb
//
// Notes:
//
// Copyright (c) 2011 R.C. Elston
// All Rights Reserved
//============================================================================

#include "assoc/Residuals.h"

namespace SAGE
{

namespace ASSOC
{

//============================================================================
// IMPLEMENTATION:  Residuals
//============================================================================
//
Residuals::Residuals()
{}

Residuals::Residuals(const vector<double>& raw_residuals, double rescaling_factor)
{
  size_t  count = raw_residuals.size();
  for(size_t i = 0; i < count; ++i)
  {
    double  residual = raw_residuals[i];
    //if(! isnan(residual))
    //{
      residual *= rescaling_factor;
    
      my_data.push_back(residual);
      my_stats.add(residual);
    //}
  }
}


}
}

