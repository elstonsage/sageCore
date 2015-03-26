#include "maxfunapi/Maximizer.h"

namespace SAGE {
namespace MAXFUN {

//=======================================================================
//  MAXFUN::maximize()
//=======================================================================
Results 
maximize(ParameterMgr& info, MaxFunction& func, const SequenceCfg& config, const DebugCfg& debug)
{
  return Maximizer::Maximize(config, info, func, debug);
}
  
//=======================================================================
//  MAXFUN::evaluateOnce()
//=======================================================================
pair<double, int>
evaluateOnce(ParameterMgr& info, MaxFunction& func, bool use_initial_ests)
{
  return Maximizer::evaluateOnce(info, func, use_initial_ests);
}

//=======================================================================
//  MAXFUN::maximizeDefault()
//=======================================================================
Results 
maximizeDefault(ParameterMgr& info, MaxFunction& func)
{
  SequenceCfg config  (SequenceCfg::DEFAULT_MAXIMIZATION);
  DebugCfg    debug   (DebugCfg::NO_DEBUG_INFO);

  return maximize(info, func, config, debug);
}   

//=======================================================================
//  MAXFUN::Maximizer::maximize()
//=======================================================================
Results
Maximizer::maximize(Function& func, const SequenceCfg& seq, const DebugCfg& debug)
{
  return Maximize(seq, func.getMgr(), func, debug);
}

//=======================================================================
//  MAXFUN::Maximizer::Maximize()
//=======================================================================
Results
Maximizer::Maximize(const SequenceCfg& seq, ParameterMgr& info, MaxFunction& func, const DebugCfg& debug)
{
  // Finalize the submodels in the ParameterMgr:
  int err = info.finalizeSubmodels();
  
  // NOTE! Must process return code!
  if(err != 0)
    SAGE_internal_error();

  // If debug is enabled, output seq info:
  if(debug.getReportBasicConfig())
    seq.dump(debug);

  // Set up local variables:
  APIMaxFunction  api_maxfunction(func, info, debug);
  Maxfun          maxfun_instance(api_maxfunction);

  // Set up maximization:
  AddParameters(info, maxfun_instance);

  // Maximize and save results:
  Results results = MaximizeFunction(info, seq, maxfun_instance, api_maxfunction, debug);
  
  // Return results object:
  return results;
}


//=======================================================================
//  MAXFUN::Maximizer::calculateScoreTestStatistic()
//=======================================================================
double
Maximizer::calculateScoreTestStatistic(Function& func, const DebugCfg& debug)
{
  ParameterMgr&  info = func.getMgr();

  // Finalize the submodels in the ParameterMgr:
  int err = info.finalizeSubmodels();
  
  // NOTE! Must process return code!
  if(err != 0)
    SAGE_internal_error();

  // Set up local variables:
  APIMaxFunction  api_maxfunction(func, info, debug);
  Maxfun          maxfun_instance(api_maxfunction);

  // Set up score test:
  AddParameters(info, maxfun_instance);
  
  cout << "\nCalling obtain score ..." << endl;

  return  maxfun_instance.obtain_score();
}


//=======================================================================
//  AddParameters()
//=======================================================================
void 
Maximizer::AddParameters(ParameterMgr& info, Maxfun& maxfun, bool use_initial_ests)
{
  // 1. Set number of parameters:

  maxfun.nt() = info.getParamCount();

  // 2. Calculate initial estimates of internal/external parameters:

  info.calculateInternalInitialEstimates();

  // 2. Add parameters:

  for(ParameterIterator p = info.getParamBegin(); p != info.getParamEnd(); p++)
    AddParameter(*p, maxfun, use_initial_ests);
}

//=======================================================================
//  AddParameter()
//=======================================================================
void 
Maximizer::AddParameter(Parameter& param, Maxfun& maxfun, bool use_initial_est)
{
  int idx = param.getIndex();

  maxfun.set_label (idx, paramname2str(param.getGroupName(), param.getName()));

  maxfun.thin      (idx) = use_initial_est ? param.getInitialEstimate() : param.getCurrentEstimate();
  maxfun.thl       (idx) = param.getLowerBound      ();
  maxfun.thu       (idx) = param.getUpperBound      ();
  maxfun.istin     (idx) = param.getInitialType     ();
  maxfun.stpin     (idx) = param.getInitialStepsizeFactor();

     if(param.scoretest) // due to JA for score test
       {  // this is a score test parameter
          maxfun.set_score_calc(idx,1); // do the score test
          param.setInitialType(Parameter::FIXED); // declaring the parameter fixed
          maxfun.istin (idx) = param.getInitialType(); // overwriting the initial type
       }
     else
       { 
          maxfun.set_score_calc(idx,0); // (skip score) due to JA for score test 
       }     
}

//=======================================================================
//  setParamIndices()
//=======================================================================
void 
Maximizer::setParamIndices(ParameterMgr& info, const Maxfun_Data& data, Results& results)
{
  // 0. Set up local variables:

  size_t param_idx = 0,
         deriv_idx = 0,
         var_idx   = 0;

  // 1. Process each parameter:

  for(ParameterIterator p = info.getParamBegin(); p != info.getParamEnd(); p++)
  {
    param_idx = p->getIndex();

    if(data.ist(param_idx) == (int)Parameter::INDEPENDENT_FUNCTIONAL ||
       data.ist(param_idx) == (int)Parameter::DEPENDENT              ||
       data.ist(param_idx) == (int)Parameter::INDEPENDENT)
      p->setDerivIndex(deriv_idx++);
    else
      p->setDerivIndex((size_t)-1);
   
    if(results.getCovMatrixStatus() == 0 && data.ist(param_idx) != (int)Parameter::FIXED)
      p->setVarIndex(var_idx++);
    else
      p->setVarIndex((size_t)-1);

    p->updateFinalOutput(data, results.getCovarianceMatrix(), results.getCovMatrixStatus());

  } // End parameter loop
}

//=======================================================================
//  MaximizeFunction()
//=======================================================================
Results
Maximizer::MaximizeFunction(ParameterMgr& info, const SequenceCfg& seq, Maxfun& maxfun, 
                                  APIMaxFunction& api_maxfunction, const DebugCfg& debug)
{
  // 0. Set up local variables:

  bool  previous_step_converged = false;
  bool  checked_initial_val     = false;
  const RunCfg*  run  = NULL;
  int  best_run = -1;

  vector<boost::shared_ptr<Results> >  results(seq.getRunCfgs().size());

  // 1. Loop through each RunCfg, transfer its contents to Maxfun, and maximize:

  for(size_t i = 0; i < seq.getRunCfgs().size(); i++)
  {

    // 1.1. Fetch RunCfg:

    run = &seq.getRunCfgs()[i];

    // 1.2. Debug header:

    if(debug.reportEachRunCfg())
    {
      std::ostringstream t;
      
      t << "MAXIMIZATION " << i << " (0 - " << seq.getRunCfgs().size() - 1 << ") " << seq.getSequenceName();
      
      debug.getOutputStream() << OUTPUT::Section(t.str());
    }
    
    // 1.3. Input debug:

    if(debug.getReportNonParametricInput())
    {
      run->Dump(debug);
    }

    // 1.3b Input parameters debug:

    if(debug.getReportParametricInput())
    {
      OUTPUT::Table t("Parametric Input");
      
      t << OUTPUT::TableColumn("Parameter") << OUTPUT::TableColumn("Initial type") << OUTPUT::TableColumn("Initial estimate"); 
     
      for(ParameterIterator p = info.getParamBegin(); p != info.getParamEnd(); p++)
      {
        OUTPUT::TableRow r = OUTPUT::TableRow()
          << p->getName()
          << ParamTypeEnum2str((Parameter::ParamTypeEnum)maxfun.istin(p->getIndex()))
          << maxfun.thin(p->getIndex());
          
        t << r;
      }
      debug.getOutputStream() << t;
    }

    // 1.4. Check the previous_step_converged status in case we need to skip this step:

    if(((run->control_option == RunCfg::PREVIOUS_NONCONVERGENCE) && (previous_step_converged == true)) ||
       ((run->control_option == RunCfg::PREVIOUS_CONVERGENCE)    && (previous_step_converged == false)))
    {
      if(debug.reportEachRunCfg())
        debug.getOutputStream() << OUTPUT::NamedString("Note", "Skipping this RunCfg...") << endl;

      continue;
    }

    // 1.5. Transfer RunCfg info to maxfun:

    seq.getRunCfgs()[i].TransferContentsToMaxfun(maxfun);

    // 1.5b. Output function evaluation header if requested:

    OUTPUT::Table table("MAXFUN RUNTIME OUTPUT");
    
    if(debug.reportEachRunCfg())
      table.enableRuntimeOutput(&debug.getOutputStream());

    if(debug.getReportEvaluations() != 0)
    {
      table << OUTPUT::TableColumn("iter.") << OUTPUT::TableColumn("f(...)");

      for(ParameterIterator p = info.getParamBegin(); p != info.getParamEnd(); p++)
        table << OUTPUT::TableColumn(p->getName());
    }

    api_maxfunction.setTable(&table);

    // 1.5a Check initial value if necessary:

    if(!checked_initial_val)
    {
      table.insertRowMessage("Initial function evaluation");

      std::pair<double, int> initial_results = evaluateOnce(info, *maxfun.fun, true);
      double                 val             = initial_results.first;
      int                    err_code        = initial_results.second;

      if(err_code || !finite(val) || SAGE::isnan(val))
      {
        table.insertRowMessage("Error: No valid initial function evaluation with these parameters.");

        if(debug.reportEachRunCfg())
          debug.getOutputStream() << table;

        Results BadResults(maxfun.get_results(), debug);

        BadResults.setValidInitialValue(false);

        return BadResults;
      }
      else
      {
        if(debug.reportEachRunCfg())
          table.insertRowMessage("Initial functional evaluation valid.");
      }

      if(debug.reportEachRunCfg())
        table.insertRowMessage("Continuing with maximization:");

      checked_initial_val = true;
    }

    // 1.6. Execute maximization:

// for the sake of checking  

       maxfun.run();

    // 1.7. Grab output:

    results[i] = boost::shared_ptr<Results>(new Results(maxfun.get_results(), debug));
    
    // 1.9. Set previous_step_converged:

    previous_step_converged = results[i]->getConverged();

    // 1.10. If this converged, set the best_run:

    if(results[i]->getConverged())
      best_run = i;

    // 1.11. Set final parameter estimates and types to initial parameter estimates and types:

    for(size_t j = 0; j < (size_t)maxfun.nt(); ++j)
    {
      maxfun.thin  (j) = maxfun.param (j);

           if(maxfun.ist(j) == (int)Parameter::INDEPENDENT_FUNCTIONAL)        maxfun.istin(j) = (int)Parameter::INDEPENDENT_FUNCTIONAL;
      else if(maxfun.ist(j) == (int)Parameter::INDEPENDENT)                   maxfun.istin(j) = (int)Parameter::INDEPENDENT;
      else if(maxfun.ist(j) == (int)Parameter::DEPENDENT)                     maxfun.istin(j) = (int)Parameter::DEPENDENT;
      else if(maxfun.ist(j) == (int)Parameter::FIXED)                         maxfun.istin(j) = (int)Parameter::FIXED;
      else if(maxfun.ist(j) == (int)Parameter::IND_FUNC_FIXED_AT_BOUND)       maxfun.istin(j) = (int)Parameter::INDEPENDENT_FUNCTIONAL;
      else if(maxfun.ist(j) == (int)Parameter::IND_FIXED_AT_BOUND)            maxfun.istin(j) = (int)Parameter::INDEPENDENT;
      else if(maxfun.ist(j) == (int)Parameter::IND_FUNC_FIXED_NEAR_BOUND)     maxfun.istin(j) = (int)Parameter::INDEPENDENT_FUNCTIONAL;
      else if(maxfun.ist(j) == (int)Parameter::IND_FIXED_NEAR_BOUND)          maxfun.istin(j) = (int)Parameter::INDEPENDENT;
      else if(maxfun.ist(j) == (int)Parameter::IND_FUNC_FIXED_NOT_NEAR_BOUND) maxfun.istin(j) = (int)Parameter::INDEPENDENT_FUNCTIONAL;
      else if(maxfun.ist(j) == (int)Parameter::IND_FIXED_NOT_NEAR_BOUND)      maxfun.istin(j) = (int)Parameter::INDEPENDENT;
    }

  }     // End RunCfg's loop

  // 2. Return either the best_run Results or the last one:
  size_t result_to_return = 0;  

  if(seq.getKeepBestRun())
  {
    debug.getOutputStream() << OUTPUT::NamedString("Note", "Keeping results for maximization #" + doub2str(best_run));

    result_to_return = best_run;
  }
  else
  {
    // We should return the final results unless for some reason (prior
    // convergence) that wasn't run, in which case, we return the best
    // we've had.
    //
    // Note that this isn't perfect.  It does not detect cases where
    // the last wasn't run due to prior *non* convergence and/or cases where
    // none of the Runs in the sequence were run.
    if(results[results.size() - 1].get() != NULL)
    {
      result_to_return = results.size() - 1;
    }
    else
    {
      result_to_return = best_run;
    }
  }

  // Assign name:
  results[result_to_return]->setSequenceName(seq.getSequenceName());
  results[result_to_return]->setFunctionName(seq.getFunctionName()); 

  // Update ParameterMgr with new indices:
  setParamIndices(info, maxfun.get_results(), *results[result_to_return]);

  // Copy ParameterMgr into Results:
  results[result_to_return]->inputParameterMgr(info);


/*
      if (results[result_to_return]->compscorestat() &&  
            results[result_to_return]->getConverged() ) 
          // for the score test, do only if convergence reached
      {

        bool goodfirstderiv = maxfun.score_deriv1();

        bool goodsecderiv = false;

        if (goodfirstderiv) goodsecderiv = maxfun.score_deriv2();

         if ( (goodsecderiv) && (goodfirstderiv) )
         {
           double score = maxfun.calc_score();
           results[result_to_return]->score = score; // inserting score into results onject
         }
      }

end of trial code
*/

  // Return the select results object:
  return *results[result_to_return];
}

//=======================================================================
//  evaluateOnce()
//=======================================================================
pair<double, int>
Maximizer::evaluateOnce(Function& func, bool use_initial_ests)
{
  return evaluateOnce(func.getMgr(), func, use_initial_ests);
}

//=======================================================================
//  evaluateOnce()
//=======================================================================
pair<double, int>
Maximizer::evaluateOnce(ParameterMgr& mgr, MaxFunction& func, bool use_initial_ests)
{
  // Finalize the submodels:
  int err = mgr.finalizeSubmodels();

  // NOTE! Must process return code!
  if(err != 0)
    SAGE_internal_error();

  // Create dummy DebugCfg:
  DebugCfg debug_cfg(DebugCfg::NO_DEBUG_INFO);

  // Build APIMaxFunction:
  APIMaxFunction my_func(func, mgr, debug_cfg);

  // Build Maxfun out of that:
  Maxfun maxfun(my_func);

  // Add the parameters to Maxfun:
  AddParameters(mgr, maxfun, use_initial_ests);

  // Now carry out a single evaluation of Maxfun:

  vector<double> temp_params(maxfun.nt());

  for(size_t j = 0; j < temp_params.size(); j++) { 
    temp_params[j] = maxfun.thin(j); 
  }

  double val      = 0.0;
  int    nfe      = 0,
         err_code = 0;

  maxfun.evaluate(temp_params, val, nfe, err_code);

  return make_pair(val, err_code);
}

} // End namespace MAXFUN
} // End namespace SAGE
