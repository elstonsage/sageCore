//============================================================================
//
//  File:  AnalysisWrapper.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//============================================================================

#include "ageon/AnalysisWrapper.h"

namespace SAGE {
namespace AO {

//=========================================================================
//
//  AnalysisWrapper(...) CONSTRUCTOR
//
//=========================================================================
AnalysisWrapper::AnalysisWrapper(
  const Model&                                 mod,
  const RPED::RefMultiPedigree&                RMP, 
  const SAMPLING::PartitionedMemberDataSample& sample,
        int                                    t,
        AnalysisOutput&                        output) 
  :
  my_RMP           (RMP),
  my_sample        (sample),
  my_analysis_type (t),
  my_model         (mod),
  my_calculator    (my_model, my_sample, my_analysis_type)
{ 
  // Output runtime info:

  cout << "Performing '"
       << my_model.get_title()
       << "' with "
       << ConvertAnalysisType(my_analysis_type, mod.get_pool_class())
       << "..." 
       << endl;

  if( t == 3 )
    cout << endl;

  // Calculate initial estimates:

  my_model.calculate_initial_est(my_sample, my_analysis_type);

  // Maximize:

  MAXFUN::SequenceCfg sequence_cfg (MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, ConvertAnalysisType(my_analysis_type, mod.get_pool_class()), "ln likelihood");

  if(my_model.getDebugCfg().hasAnyReportedOutput())
    my_model.getDebugCfg().setDebugOutput(my_model.get_title() + "-" + ConvertAnalysisType(my_analysis_type, mod.get_pool_class()) + ".max", false);

  MAXFUN::Results re = MAXFUN::maximize(my_model.GetParameterMgr(), my_calculator, sequence_cfg, my_model.getDebugCfg());

#if 0
  // May need these codes in for debug = true case.
  
  MAXFUN::DebugCfg db(MAXFUN::DebugCfg::COMPLETE);

  MAXFUN::Results re = MAXFUN::Maximizer::Maximize(sequence_cfg, my_model.GetParameterMgr(), my_calculator, db);

  std::cout << MAXFUN::OutputFormatter::convertEstimates (re)
            << MAXFUN::OutputFormatter::convertMatrix    (re)
            << MAXFUN::JointTest(re, re).summarizeAsTable();
#endif

  output.input(re, my_analysis_type);
  output.input(my_model, my_analysis_type);
  output.input(my_calculator.get_mcc(), my_analysis_type);
}

} // End namespace AO
} // End namespace SAGE
