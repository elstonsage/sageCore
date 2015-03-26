#ifndef AO_CALCULATOR_H
#define AO_CALCULATOR_H
//======================================================================
//
//  File:	Calculator.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================


#include <string>
#include "maxfun/maxfun.h"
#include "sampling/sampling.h"
#include "ageon/Datatypes.h"
#include "ageon/MemberCovariateCalculator.h"
#include "ageon/Model.h"
#include "ageon/Kernel.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
//  class Calculator
//
//======================================================================
class Calculator : public MaxFunction
{
  public:
    //==================================================================
    // Constructor:
    //==================================================================

    Calculator(Model & mod, const SAMPLING::PartitionedMemberDataSample & sample, int t);

    //==================================================================
    // Public utility functions:
    //==================================================================

    virtual double evaluate      (vector<double> & params);
    virtual int    update_bounds (vector<double> & params);

    //==================================================================
    // Public accessors:
    //==================================================================

    const MemberCovariateCalculator & get_mcc() const;

  private:
    //==================================================================
    // Data members:
    //==================================================================

    Model      & my_model;
    Kernel       my_kernel;
    int          my_analysis_type;
};

} // End namespace AO
} // End namespace SAGE

#endif
