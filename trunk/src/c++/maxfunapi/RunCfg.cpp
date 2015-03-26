//=======================================================================
//
//  File:	RunCfg.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "globals/SAGEConstants.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/RunCfg.h"
#include "LSF/parse_ops.h"

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
// RunCfg()
//=======================================================================
RunCfg::RunCfg(MaximizationMethodEnum _method, int _max_iterations)
{
  // 1. Assign required values:

	method         = _method;
	max_iterations = _max_iterations;

  // 2. Set all others to unspecified:

	epsilon1       = QNAN;
	epsilon2       = QNAN;
	epsilon_delta  = QNAN;
	epsilon_trunc  = QNAN;
	stepsize       = QNAN;

	second_deriv   = SEC_DERIV_UNSPECIFIED;
	var_cov        = VAR_COV_UNSPECIFIED;

  // 3. Set control feature:

	control_option = UNCONDITIONAL;
}

//=======================================================================
// RunCfg()
//=======================================================================
RunCfg::RunCfg(const RunCfg & other)
{
  copy(other);
}

//=======================================================================
// operator=()
//=======================================================================
RunCfg& 
RunCfg::operator=(const RunCfg& other)
{
  if(&other != this)
  {
    copy(other);
  }

  return *this;
}

//=======================================================================
// copy()
//=======================================================================
void
RunCfg::copy(const RunCfg& other)
{
  control_option = other.control_option;
  epsilon1       = other.epsilon1;
  epsilon2       = other.epsilon2;
  epsilon_delta  = other.epsilon_delta;
  epsilon_trunc  = other.epsilon_trunc;
  max_iterations = other.max_iterations;
  method         = other.method;
  second_deriv   = other.second_deriv;
  stepsize       = other.stepsize;
  var_cov        = other.var_cov;
}

//=======================================================================
// TransferContentsToMaxfun()
//=======================================================================
void 
RunCfg::TransferContentsToMaxfun(Maxfun& maxfun) const
{
  // 1. Required components:

	maxfun.method () = method;
	maxfun.maxit  () = maxfun.nt() * max_iterations;

  // 2. Optional components:

	if(!SAGE::isnan(epsilon1))      maxfun.epsc1 () = epsilon1;
	if(!SAGE::isnan(epsilon2))      maxfun.epsc2 () = epsilon2;
	if(!SAGE::isnan(epsilon_delta)) maxfun.epsd  () = epsilon_delta;
	if(!SAGE::isnan(epsilon_trunc)) maxfun.epst  () = epsilon_trunc;
	if(!SAGE::isnan(second_deriv))  maxfun.ihit  () = second_deriv;
	if(!SAGE::isnan(stepsize))      maxfun.yota  () = stepsize;
	if(!SAGE::isnan(var_cov))       maxfun.ixvc  () = var_cov;
}

//=======================================================================
// Dump()
//=======================================================================
void
RunCfg::Dump(const DebugCfg& debug) const
{
  // Header

  OUTPUT::Table t("MAXFUN INPUT");

  // Method

  OUTPUT::TableRow method_row;
  
  method_row << "method";
  
  switch(method)
  {
    case COMPLETE_DIRECT       : method_row << "Complete direct search";                                                                 break;
    case DIRECT_WITHOUT        : method_row << "Direct search without 2^Ni-trial search";                                                break;
    case NEWTON_WITHOUT_SECOND : method_row << "Newton-Raphson method, without repeated calculation of second derivatives";              break;
    case NEWTON_WITH_SECOND    : method_row << "Newton-Raphson method with second derivatives recalculated for each iteration";          break;
    case VAR_METRIC_IDENTITY   : method_row << "Variable metric method, with initial B = identity";                                      break;
    case VAR_METRIC_ESTIMATE   : method_row << "Variable metric method, with initial B = -H calculated at initial parameter estimates."; break;
  }

  t << method_row;

  // Control option

  OUTPUT::TableRow control_opt_row;
  
  control_opt_row << "control_option";

  switch(control_option)
  {
    case UNCONDITIONAL           : control_opt_row << "unconditional, Do unconditionally."; break;
    case PREVIOUS_NONCONVERGENCE : control_opt_row << "previous_nonconvergence, Do only if previous procedure didn't reach convergence."; break;
    case PREVIOUS_CONVERGENCE    : control_opt_row << "previous_convergence, Do only if previous procedure did reach converge."; break;
  }

  t << control_opt_row;
  
  t << (OUTPUT::TableRow() << "max_iterations" << (!SAGE::isnan(max_iterations) ? doub2str(max_iterations) : "Must be specified"));
  
  t << OUTPUT::NamedString("Note", "max_iterations is actually multiplied by the total number of parameters before maximization begins.");

  t << (OUTPUT::TableRow() << "epsilon1"      << (SAGE::isnan(epsilon1)      ? "Default (see documentation)" : OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::SCIENTIFIC).render(OUTPUT::Double(epsilon1))))
    << (OUTPUT::TableRow() << "epsilon2"      << (SAGE::isnan(epsilon2)      ? "Default (see documentation)" : OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::SCIENTIFIC).render(OUTPUT::Double(epsilon2))))
    << (OUTPUT::TableRow() << "epsilon_delta" << (SAGE::isnan(epsilon_delta) ? "Default (see documentation)" : OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::SCIENTIFIC).render(OUTPUT::Double(epsilon_delta))))
    << (OUTPUT::TableRow() << "epsilon_trunc" << (SAGE::isnan(epsilon_trunc) ? "Default (see documentation)" : OUTPUT::RenderingRules(OUTPUT::HasNumberFormat::SCIENTIFIC).render(OUTPUT::Double(epsilon_trunc))));

  // Var-cov option

  OUTPUT::TableRow var_cov_row;
  
  var_cov_row << "var_cov";

  switch(var_cov)
  {
    case VAR_COV_UNSPECIFIED : var_cov_row << "Default (see documentation)"; break;
    case NO_MATRIX           : var_cov_row << "Not computed beyond what is necessary for parameter estimation."; break;
    case INITIAL             : var_cov_row << "Computed for final estimates unless one is available that is no more than one iteration old."; break;
    case FINAL               : var_cov_row << "Computed for final estimates (using iterative method)."; break;
  }

  t << var_cov_row;

  // Sec-deriv option

  OUTPUT::TableRow deriv_row;
  
  deriv_row << "second_deriv";

  switch(second_deriv)
  {
    case SEC_DERIV_UNSPECIFIED : deriv_row << "Default (see documentation)"; break;
    case SINGLE                : deriv_row << "Use a single approximation"; break;
    case ITERATIVE             : deriv_row << "Estimate through an iterative process"; break;
  }

  t << deriv_row;
 
  // Stepsize

  t << (OUTPUT::TableRow() << "stepsize" << (!SAGE::isnan(stepsize) ? doub2str(stepsize) : "Default (see documentation)"));

  debug.getOutputStream() << t;
}

}
} 
