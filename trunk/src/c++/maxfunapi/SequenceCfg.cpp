//=======================================================================
//
//  File:	SequenceCfg.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "maxfunapi/SequenceCfg.h"

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
// CONSTRUCTOR
//=======================================================================
SequenceCfg::SequenceCfg(SequenceTemplateEnum sequence_template, string sequence_name, string function_name)
{
  my_sequence_name = sequence_name;
  my_function_name = function_name;
  my_KeepBestRun   = false;
  my_RunCfgs.clear();

  setSequenceTemplate(sequence_template);
}

//=======================================================================
// COPY CONSTRUCTOR
//=======================================================================
SequenceCfg::SequenceCfg(const SequenceCfg& other)
{
  copy(other);
}

//=======================================================================
// operator=()
//=======================================================================
SequenceCfg& 
SequenceCfg::operator=(const SequenceCfg& other)
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
SequenceCfg::copy(const SequenceCfg& other)
{
  my_KeepBestRun       = other.my_KeepBestRun;
  my_RunCfgs           = other.my_RunCfgs;
  my_sequence_template = other.my_sequence_template;
  my_sequence_name     = other.my_sequence_name;
  my_function_name     = other.my_function_name;
}

//=======================================================================
// setKeepBestRun()
//=======================================================================
int 
SequenceCfg::setKeepBestRun(bool keep)
{
  my_KeepBestRun = keep; return 0;
}

//=======================================================================
// getKeepBestRun()
//=======================================================================
bool
SequenceCfg::getKeepBestRun() const
{
  return my_KeepBestRun;
}

//=======================================================================
// setSequenceTemplate(SequenceTemplateEnum sequence_template)
//=======================================================================
void
SequenceCfg::setSequenceTemplate(SequenceTemplateEnum sequence_template)
{
  // 1. Assign the SequenceTemplateEnum:

	my_sequence_template = sequence_template;

  // 2. Process the sequence template:

	switch(my_sequence_template)
	{

  // 2.1. User-defined:

	  case SequenceCfg::USER_DEFINED: 
	                    break;

  // 2.2. Default:

	  case SequenceCfg::DEFAULT_MAXIMIZATION: 

                            addRunCfg(RunCfg::DIRECT_WITHOUT, 1);
                            getLatestRunCfg().epsilon1 = 1e-3;
                            getLatestRunCfg().epsilon2 = 1e-10;

                            addRunCfg(RunCfg::VAR_METRIC_IDENTITY, 20);
                            getLatestRunCfg().epsilon1       = 1e-3;  
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;
                            getLatestRunCfg().var_cov        = RunCfg::FINAL;

                            addRunCfg(RunCfg::DIRECT_WITHOUT, 50);
                            getLatestRunCfg().epsilon1       = 1e-4;
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;

                            addRunCfg(RunCfg::VAR_METRIC_ESTIMATE, 20);
                            getLatestRunCfg().epsilon1       = 1e-10;
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().var_cov        = RunCfg::FINAL;

                            addRunCfg(RunCfg::VAR_METRIC_ESTIMATE, 20);
                            getLatestRunCfg().epsilon1       = 1e-8;
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().var_cov        = RunCfg::FINAL;
                            getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;

                            addRunCfg(RunCfg::VAR_METRIC_ESTIMATE, 20);
                            getLatestRunCfg().epsilon1       = 1e-6;
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().var_cov        = RunCfg::FINAL;
                            getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;

                            addRunCfg(RunCfg::VAR_METRIC_ESTIMATE, 20);
                            getLatestRunCfg().epsilon1       = 1e-4;
                            getLatestRunCfg().epsilon2       = 1e-10;
                            getLatestRunCfg().var_cov        = RunCfg::FINAL;
                            getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;

                            break;
	}
}

//=======================================================================
// getSequenceTemplate()
//=======================================================================
SequenceCfg::SequenceTemplateEnum
SequenceCfg::getSequenceTemplate()
{
  return my_sequence_template;
}

//=======================================================================
// getSequenceName()
//=======================================================================
const std::string&
SequenceCfg::getSequenceName() const
{
  return my_sequence_name;
}

//=======================================================================
// getFunctionName()
//=======================================================================
const std::string&
SequenceCfg::getFunctionName() const
{
  return my_function_name;
}

void
SequenceCfg::setFunctionName(const std::string & name)
{
  my_function_name = name;
}

void
SequenceCfg::setSequenceName(const std::string& name)
{
  my_sequence_name = name;
}

//=======================================================================
// Dump()
//=======================================================================
void
SequenceCfg::dump(const DebugCfg& debug) const
{
  // 1. Maxtype:

	debug.getOutputStream() << "SequenceTemplateEnum = " << my_sequence_template << ", ";

	switch(my_sequence_template)
	{
	  case USER_DEFINED         : debug.getOutputStream() << "User-defined" << endl; break;
	  case DEFAULT_MAXIMIZATION : debug.getOutputStream() << "Default"      << endl; break;
	}

	if(my_KeepBestRun)
	  debug.getOutputStream() << "Maxfun will keep the best run." << endl;

	debug.getOutputStream() << endl;
}

//=======================================================================
// AddRunCfg()
//=======================================================================
int 
SequenceCfg::addRunCfg(RunCfg::MaximizationMethodEnum method, int max_iterations)
{
  my_RunCfgs.push_back(RunCfg(method, max_iterations));

  return 0;
}

//=======================================================================
// DuplicateLatestRunCfg()
//=======================================================================
int 
SequenceCfg::duplicateLatestRunCfg()
{
  int i = my_RunCfgs.size() - 1;

  my_RunCfgs.push_back(my_RunCfgs[i]);

  return 0;
}

//=======================================================================
// GetLatestRunCfg()
//=======================================================================
RunCfg & 
SequenceCfg::getLatestRunCfg()
{
  if(my_RunCfgs.size())
  {
    return my_RunCfgs[my_RunCfgs.size() - 1];
  }
  else
  {
    addRunCfg(RunCfg::DIRECT_WITHOUT, 1);

    return getLatestRunCfg();
  }
}

//=======================================================================
// GetRunCfg()
//=======================================================================
RunCfg& 
SequenceCfg::getRunCfg(int i)
{
  return my_RunCfgs[i];
}

//=======================================================================
// GetRunCfgs()
//=======================================================================
const vector<RunCfg> & SequenceCfg::getRunCfgs() const { return my_RunCfgs; }
      vector<RunCfg> & SequenceCfg::getRunCfgs()       { return my_RunCfgs; }

}
} 
