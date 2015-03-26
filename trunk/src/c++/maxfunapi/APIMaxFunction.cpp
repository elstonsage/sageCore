#include "maxfunapi/APIMaxFunction.h"

namespace SAGE   {
namespace MAXFUN {

//======================================================================
// CONSTRUCTOR
//======================================================================
APIMaxFunction::APIMaxFunction(SAGE::MaxFunction & func, ParameterMgr & info, const DebugCfg & debug) :
	my_max_function   (func),
	my_parameter_mgr  (info),
	my_debug_cfg      (debug)
{
  nfe      = 0;
  my_table = NULL;
}

//======================================================================
// COPY CONSTRUCTOR
//======================================================================
APIMaxFunction::APIMaxFunction(const APIMaxFunction & other) :
        SAGE::MaxFunction (other),
	my_max_function   (other.my_max_function),
	my_parameter_mgr  (other.my_parameter_mgr),
	my_debug_cfg      (other.my_debug_cfg)
{
  nfe      = other.nfe;
  my_table = other.my_table;
}

//============================================
//  setTable(...)
//============================================
void 
APIMaxFunction::setTable(OUTPUT::Table * table)
{
  my_table = table;
}
    

//======================================================================
// evaluate(...)
//======================================================================
double 
APIMaxFunction::evaluate(vector<double> & params)
{ 
     if( (my_debug_cfg.getType() == SAGE::MAXFUN::DebugCfg::COMPLETE) && (my_debug_cfg.iter_end) )
       { 
         if (Maxfun::iteration_ended) output_iteration_end_results();
       }

  // 1. Increment number of function evaluations:
        
	nfe++;

  // 2. Fetch result:

	double val = my_max_function(params);
        lastvalue = val;

  // 3. Output result if requested: 

	if((my_debug_cfg.getReportEvaluations() != 0) && ((nfe % my_debug_cfg.getReportEvaluations()) == 0) \
            && (!(my_debug_cfg.iter_end)) )
	{
          OUTPUT::TableRow r;
	  
	  r << nfe << val;

	  for(ParameterIterator p = getParameterMgr().getParamBegin(); p != getParameterMgr().getParamEnd(); p++)
            r << p->getCurrentEstimate();
          
          if(my_table != NULL)
            (*my_table) << r;
	}

  // 4. Return result:

        return val;
}

//======================================================================
// update_bounds(...)
//======================================================================
int
APIMaxFunction::update_bounds(vector<double> & params)
{
  // Pre-precess - Initialize error code

        int err_code = 0;

  // 0. Update ParameterMgr object:

	err_code = getParameterMgr().update(params);

	if(err_code)
	  return err_code;

  // 1. Invoke update_bounds on APIMaxFunction:

	err_code = my_max_function.depar(params);

	if(err_code)
	  return err_code;

  // 2. Update dependents:

	err_code = getParameterMgr().updateDependents(params);

	if(err_code)
	  return err_code;

  // 3. Return success:

	return 0;
}

//======================================================================
// getParameterMgr()
//======================================================================
ParameterMgr & 
APIMaxFunction::getParameterMgr() const
{
  return my_parameter_mgr;
}

void APIMaxFunction::output_iteration_end_results()
{ // so we print out only at the end of a maxfun iteration
	if(my_debug_cfg.getReportEvaluations() != 0) 
	{
	  OUTPUT::TableRow r;
	  
          r << nfe << lastvalue;

	  for(ParameterIterator p = getParameterMgr().getParamBegin(); p != getParameterMgr().getParamEnd(); p++)
            r << p->getCurrentEstimate();
          
          if(my_table != NULL)
            (*my_table) << r;

        }

       Maxfun::iteration_ended = false;

       return;
}
}} // End namespace

