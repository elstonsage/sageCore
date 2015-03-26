#ifndef AO_ANALYSIS_WRAPPER_H
#define AO_ANALYSIS_WRAPPER_H
//=============================================================================
//
//  File:	AnalysisWrapper.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=============================================================================


#include <fstream>
#include "LSF/parse_ops.h"
#include "app/output_streams.h"
#include "rped/rped.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/maxfunapi.h"
#include "sampling/sampling.h"
#include "ageon/AnalysisOutput.h"
#include "ageon/Datatypes.h"
#include "ageon/Model.h"
#include "ageon/Calculator.h"

namespace SAGE {
namespace AO   {

//=============================================================================
//
// class AnalysisWrapper
//
//=============================================================================
class AnalysisWrapper
{
  public:
    //=========================================================================
    // Constructor:
    //=========================================================================
  
    AnalysisWrapper(
      const Model &, 
      const RPED::RefMultiPedigree &, 
      const SAMPLING::PartitionedMemberDataSample & sample, 
      int t, 
      AnalysisOutput &);
    
  protected:
    //=========================================================================
    // Private mutators:
    //=========================================================================

    MAXFUN::ParameterMgr & GetParameterMgr();

    //=========================================================================
    // Private utility functions:
    //=========================================================================

    void build_calculator     ();
    void setup_maximization   ();
    void check_initial_lh     (vector<double> &);
    void perform_maximization ();
    void adjust_model         ();
    void update_output        ();
    void cleanup              ();

    //=========================================================================
    // Data members:
    //=========================================================================

    const RPED::RefMultiPedigree          & my_RMP;
    const SAMPLING::PartitionedMemberDataSample & my_sample;
	  int                                     my_analysis_type;
          Model                                   my_model;
          Calculator                              my_calculator;
};

inline MAXFUN::ParameterMgr     & AnalysisWrapper::GetParameterMgr() { return my_model.GetParameterMgr(); }

}} // End namespace

#endif
