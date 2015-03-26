//=======================================================================
//
//  File:	ParameterMgr.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//  ParameterMgr() CONSTRUCTOR
//=======================================================================
ParameterMgr::ParameterMgr()
{
  // Reset everything:
  reset();
}

//=======================================================================
//  ParameterMgr() COPY CONSTRUCTOR
//=======================================================================
ParameterMgr::ParameterMgr(const ParameterMgr & other)
{
  copy(other);
}

//=======================================================================
//  operator=()
//=======================================================================
ParameterMgr& 
ParameterMgr::operator=(const ParameterMgr& other)
{
  if(&other != this)     // Added 7-19-7. djb
  {
    copy(other);
  }

  return *this;
}

//=======================================================================
//  DESTRUCTOR
//=======================================================================
ParameterMgr::~ParameterMgr()
{
  // Make sure all submodels are informed that we're done with them.

  removeSubmodels();
}

//=======================================================================
//  copy()
//=======================================================================
void
ParameterMgr::copy(const ParameterMgr& other)
{
  // We may not copy when either submodel is in use.

  if(other.isInUse())
    SAGE_internal_error();

  if(isInUse())
    SAGE_internal_error();

  // Copy all the data structures

  group_ordering       = other.group_ordering;
  groups               = other.groups;
  param_lookup_table   = other.param_lookup_table;
  params               = other.params;
  my_transformer_infos = other.my_transformer_infos;
  my_param_calcs       = other.my_param_calcs;
  
  my_new_submodels.clear();

  for(size_t i = 0; i < other.my_new_submodels.size(); ++i)
  {
    NewSubmodelShPtr sm = other.my_new_submodels[i]->clone();
    sm->setParameterMgr(this);
    my_new_submodels.push_back(sm);
  }

}

//=======================================================================
//  copyInUseParameterMgr()
//=======================================================================
void
ParameterMgr::copyInUseParameterMgr(const ParameterMgr& other)
{
  // Copy all the data structures

  group_ordering       = other.group_ordering;
  groups               = other.groups;
  param_lookup_table   = other.param_lookup_table;
  params               = other.params;
  my_transformer_infos = other.my_transformer_infos;
  my_param_calcs       = other.my_param_calcs;

  my_new_submodels.clear();

  for(size_t i = 0; i < other.my_new_submodels.size(); ++i)
  {
    NewSubmodelShPtr sm = other.my_new_submodels[i]->clone();
    sm->setParameterMgr(this);
    my_new_submodels.push_back(sm);
  }
}

//=======================================================================
//  reset()
//=======================================================================
void
ParameterMgr::reset()
{
  // 0. We may not Reset when the ParameterMgr is actively in use.

	if(isInUse())
	  SAGE_internal_error();

  // 1. Clear data:

  my_param_calcs.clear();
	my_submodels.clear();
  my_new_submodels.clear();
	params.clear();
	param_lookup_table.clear();
	group_ordering.clear();
	groups.clear();
  my_transformer_infos.clear();

  // 2. Add the "ALL" group:

	addGroup("ALL");
}

//=======================================================================
//  operator()() #1 NON CONST
//=======================================================================
double&
ParameterMgr::operator() (string group_name, string param_name)
{
  return getParameter(group_name, param_name).operator()();
}

//=======================================================================
//  operator()() #1 CONST
//=======================================================================
double
ParameterMgr::operator() (string group_name, string param_name) const
{
  return getParameter(group_name, param_name).getCurrentEstimate();
}

//=======================================================================
//  operator()() #2 NON CONST
//=======================================================================
double&
ParameterMgr::operator() (int param_id)
{
  return params[param_id].operator()();
}

//=======================================================================
//  operator()() #2 CONST
//=======================================================================
double
ParameterMgr::operator() (int param_id) const
{
  return params[param_id].getCurrentEstimate();
}

//=======================================================================
//  getEst()
//=======================================================================
double
ParameterMgr::getEst(int param_id) const
{
  return params[param_id].getCurrentEstimate();
}

//=======================================================================
//  getParameter() #1 NON CONST
//=======================================================================
Parameter& 
ParameterMgr::getParameter(string group_name, string param_name)
{
  return params[getParamID(group_name, param_name)];
}

//=======================================================================
//  getParameter() #1 CONST
//=======================================================================
const Parameter& 
ParameterMgr::getParameter(string group_name, string param_name) const
{
  return params[getParamID(group_name, param_name)];
}

//=======================================================================
//  getParameter() #2 NON-CONST
//=======================================================================
Parameter & 
ParameterMgr::getParameter(int param_id)
{
  return params[param_id];
}

//=======================================================================
//  getParameter() #2 CONST
//=======================================================================
const Parameter& 
ParameterMgr::getParameter(int param_id) const
{
  return params[param_id];
}

//=======================================================================
//  getParameter() #3 NON-CONST
//=======================================================================
Parameter& 
ParameterMgr::getParameter(string group_name, int group_id)
{
  return params[getParamID(group_name, group_id)];
}

//=======================================================================
//  getParameter() #3 CONST
//=======================================================================
const Parameter& 
ParameterMgr::getParameter(string group_name, int group_id) const
{
  return params[getParamID(group_name, group_id)];
}

//=======================================================================
//  doesGroupExist()
//=======================================================================
bool 
ParameterMgr::doesGroupExist(string group_name) const
{
  return groups.find(group_name) != groups.end();
}

//=======================================================================
//  getGroup() NON CONST
//=======================================================================
param_group& 
ParameterMgr::getGroup(string group_name)
{
  if(doesGroupExist(group_name))
  {
    return groups.find(group_name)->second;
  }
  else
  {
    addGroup(group_name);

    return getGroup(group_name);
  }
}

//=======================================================================
//  getGroup() CONST
//=======================================================================
const param_group& 
ParameterMgr::getGroup(string group_name) const
{
  if(doesGroupExist(group_name))
  {
    return groups.find(group_name)->second;
  }
  else
  {
    addGroup(group_name);

    return getGroup(group_name);
  }
}

//=======================================================================
//  addGroup()
//=======================================================================
void
ParameterMgr::addGroup(string group_name) const
{
  map<string, param_group>::iterator g = groups.find(group_name);

  if(g == groups.end())
  {
    groups.insert(make_pair(group_name, param_group(0)));
    group_ordering.push_back(group_name);
  }
}

//=======================================================================
//  doesParamExist()
//=======================================================================
bool 
ParameterMgr::doesParamExist(string group_name, string param_name) const
{
  map<pair<string, string>, int>::const_iterator i = param_lookup_table.find(make_pair(group_name, param_name));
        
  if(i == param_lookup_table.end())
    return false;
  else
    return true;
}

//=======================================================================
//  getOrderedGroupNames()
//=======================================================================
vector<string> 
ParameterMgr::getOrderedGroupNames() const
{
  return group_ordering;
}

//=======================================================================
//  addParameter() #2
//=======================================================================
int 
ParameterMgr::addParameter(ParameterInput& param)
{
  param.index = addParameter(param.group_name,   param.param_name, 
                             param.initial_type, param.initial_estimate, 
                             param.lower_bound,  param.upper_bound);

  return param.index;
}

Parameter&
ParameterMgr::addParameterAlt(ParameterInput& param)
{
  addParameter(param);
  
  return  params.back();
}

//=======================================================================
//  addParameter() #1
//=======================================================================
int 
ParameterMgr::addParameter(string                   group_name,
                           string                   param_name,
                           Parameter::ParamTypeEnum initial_type,
                           double                   initial_estimate,
                           double                   lower_bound,
                           double                   upper_bound,
                           double                   initial_stepsize)
{
  // 0. Set up local variables:

	Parameter param;

  // 0b. Check the group_name is valid:

	if(group_name == "ALL")
	{
	  cout << "Error: Cannot add " 
	       << param_name 
	       << " with group_name = " 
	       << group_name 
	       << endl
	       << "'ALL' is a reserved group name." 
	       << endl;

	  exit(0);
	}

  // 0c. Add group if necessary:

	addGroup(group_name);

  // 1. Create parameter:

	param.setName            (param_name);
	param.setNameAbbr        (param_name);
	param.setGroupName       (group_name);
	param.setGroupIndex      (getParamCount(group_name));
	param.setInitialEstimate (initial_estimate);
	param.setInitialType     (initial_type);
	param.setLowerBound      (lower_bound);
	param.setUpperBound      (upper_bound);
	param.setIndex           (getParamCount());
	param.setInitialStepsizeFactor(initial_stepsize);

  // 2. Stick it on to the parameter vector:

	params.push_back(param);

  // 3. Stick it on to the "ALL" and specific group vectors:

	getGroup("ALL").push_back(param.getIndex());
	getGroup(group_name).push_back(param.getIndex());

  // 3. Add lookup by name:

	param_lookup_table.insert(make_pair(make_pair(group_name, param_name), param.getIndex()));

  // 4. Return index number:
  
	return param.getIndex();
}

Parameter& 
ParameterMgr::addParameterAlt(string                   group_name,
                              string                   param_name,
                              Parameter::ParamTypeEnum initial_type,
                              double                   initial_estimate,
                              double                   lower_bound,
                              double                   upper_bound,
                              double                   initial_stepsize)
{
  addParameter(group_name, param_name, initial_type, 
               initial_estimate, lower_bound, upper_bound, initial_stepsize);
               
  return  params.back();
}

//======================================================================
//  addParamCalculator()
//======================================================================
int 
ParameterMgr::addParamCalculator(
  int idx,
  ParamCalculatorShPtr calculator)
{
  my_param_calcs.insert(make_pair(idx, calculator));

  return 0;
}

//=======================================================================
//  addTransformer()
//=======================================================================
int
ParameterMgr::addTransformer(
  TransformerShCstPtr transformer, 
  const string & maximization_group_name, 
  const string & maximization_param_name,
  const string & reported_group_name, 
  const string & reported_param_name)
{
  // Ensure that the parameter in question exists:
  if(doesParamExist(reported_group_name,     reported_param_name)     == false ||
     doesParamExist(maximization_group_name, maximization_param_name) == false)
  {
    return 1;
  }

  // Get the parameters and set them up:
  Parameter & reported_param     = getParameter(reported_group_name,     reported_param_name);
  Parameter & maximization_param = getParameter(maximization_group_name, maximization_param_name);

  reported_param     .setInitialType (Parameter::DEPENDENT);
  maximization_param .setInitialType (Parameter::INDEPENDENT_FUNCTIONAL);
  maximization_param .setLowerBound  (transformer->getLowerBound());
  maximization_param .setUpperBound  (transformer->getUpperBound());

  // Add the transformer entry into the internal vector:
  TransformerInfo info;

  info.reported_param_id     = reported_param     .getIndex();
  info.maximization_param_id = maximization_param .getIndex();
  info.transformer           = transformer;

  my_transformer_infos.push_back(info);

  // Return success:
  return 0;
}

//=======================================================================
//  addSubModel()
//=======================================================================
int 
ParameterMgr::addSubModel(Submodel* sub_mod)
{
  // 0. Set up local variables:

	int err = 1;

  // 1. Check the submodel for validity.  These tests should never fail (thus
  //    they fail spectacularly)

	if(!sub_mod)
	  SAGE_internal_error();

  // 2. Check submodel for linkage:

	if(sub_mod->isLinked())
	  SAGE_internal_error();

  // 3. Link the submodel:

	sub_mod->linkToParameterMgr(this);

  // 4. Finalize submodel configuration:

	err = sub_mod->finalizeConfiguration();

	if(err) return err;

  // 5. Now add the parameters:

	err = sub_mod->addParametersToParameterMgr();

	if(err) return err;

  // 6. Set any necessary advanced parameter options:

	sub_mod->setAdvancedParameterOptions();

	my_submodels.push_back(sub_mod);

  // 7. Return success:

	return 0;
}

//=======================================================================
//  finalizeSubmodels()
//=======================================================================
int
ParameterMgr::finalizeSubmodels()
{
  int err = 0;

  for(size_t i = 0; i < my_new_submodels.size(); ++i)
  {
    if(my_new_submodels[i]->isFinalized())
      continue;

    err = my_new_submodels[i]->finalizeConfiguration();

    if(err) return err;

    err = my_new_submodels[i]->addParametersToMgr();

    if(err) return err;
    
    my_new_submodels[i]->setAdvancedParameterOptions();

    my_new_submodels[i]->setFinalized(true);
  }

  return 0;
}



//=======================================================================
//  removeSubModel()
//=======================================================================
void 
ParameterMgr::removeSubmodel(Submodel* sub_mod)
{
  // Find the submodel in the list and remove it.

  for(list<Submodel*>::iterator i = my_submodels.begin(); i != my_submodels.end(); ++i)
  {
    if(*i == sub_mod)
    {
      my_submodels.erase(i);
      return;
    }
  }

  // If we reach this point, we didn't find the submodel in the list.  This
  // should never happen, so we exit spectacularly.

  SAGE_internal_error();
}

//=======================================================================
//  removeSubModels()
//=======================================================================
void 
ParameterMgr::removeSubmodels()
{
  for(list<Submodel*>::iterator i = my_submodels.begin(); i != my_submodels.end(); ++i)
    (*i)->unlinkFromParameterMgr();

  my_submodels.clear();
}

//=======================================================================
//  update()
//=======================================================================
int 
ParameterMgr::update(vector<double>& params)
{
  // 1. Copy values out of maxfun:

	for(size_t i = 0; i < this->params.size(); i++)
	  this->params[i].setCurrentEstimate(params[i]);

  // 2. Do the transformations:

        transformParameters();

  // 1b. Update submodels:

	for(list<Submodel *>::iterator i = my_submodels.begin(); i != my_submodels.end(); ++i)
        {
          int ret_val = (*i)->update();

          if(ret_val)
            return ret_val;
        }

  // 1b. Update NEW submodels:

        for(size_t i = 0; i < my_new_submodels.size(); ++i)
        {
          int ret_val = my_new_submodels[i]->update();

          if(ret_val)
            return ret_val;
        }

  // 2. Return success:

	return 0;
}

//=======================================================================
//  transformParameters()
//=======================================================================
void
ParameterMgr::transformParameters()
{
  for(vector<TransformerInfo>::iterator transformer_itr  = my_transformer_infos.begin ();
                                        transformer_itr != my_transformer_infos.end   (); ++transformer_itr)
  {
    Parameter & reported_param     = params[transformer_itr->reported_param_id],
              & maximization_param = params[transformer_itr->maximization_param_id];

//    if(reported_param.getInitialType () != Parameter::FIXED &&
//       reported_param.isInUse        () == true)
    if(reported_param.getInitialType () != Parameter::FIXED)
    {
      double current_maximization_est = maximization_param.getCurrentEstimate(),
             new_reported_est         = transformer_itr->transformer->transformToReportedScale(this, current_maximization_est);

      reported_param.setCurrentEstimate(new_reported_est);
    }
  }
}

//=======================================================================
//  calculateInternalInitialEstimates()
//=======================================================================
void 
ParameterMgr::calculateInternalInitialEstimates()
{
  for(vector<TransformerInfo>::iterator transformer_itr  = my_transformer_infos.begin ();
                                        transformer_itr != my_transformer_infos.end   (); ++transformer_itr)
  {
    // Grab the parameters in question:

    Parameter & reported_param     = params[transformer_itr->reported_param_id],
              & maximization_param = params[transformer_itr->maximization_param_id];

    // If the reported param is FIXED, or not IN USE, then take the maximization
    // param out of the game!

//    if(reported_param.getInitialType () == Parameter::FIXED ||
//       reported_param.isInUse        () == false)
    if(reported_param.getInitialType () == Parameter::FIXED)
    {
      continue; // Added
//      maximization_param.setInUse(false);
    }

    // Otherwise, set up the initial estimate and types of the
    // maximization/reported params:

    else if(reported_param.getInitialType() == Parameter::INDEPENDENT ||
            reported_param.getInitialType() == Parameter::INDEPENDENT_FUNCTIONAL)
    {
      if(SAGE::isnan(maximization_param.getInitialEstimate()))
      {
        double reported_initial_est         = reported_param.getInitialEstimate(),
               new_maximization_initial_est = transformer_itr->transformer->transformToMaximizationScale(this, reported_initial_est);

        maximization_param.setInitialEstimate (new_maximization_initial_est);
      }
      else if(SAGE::isnan(reported_param.getInitialEstimate()))
      {
        double maximization_initial_est = maximization_param.getInitialEstimate(),
               new_reported_initial_est = transformer_itr->transformer->transformToReportedScale(this, maximization_initial_est);

        reported_param.setInitialEstimate(new_reported_initial_est);
      }

      maximization_param .setInitialType(Parameter::INDEPENDENT_FUNCTIONAL);
      reported_param     .setInitialType(Parameter::DEPENDENT);
    }
  }
}

//=======================================================================
//  updateDependents()
//=======================================================================
int 
ParameterMgr::updateDependents(vector<double>& params)
{
  for(size_t i = 0; i < this->params.size(); i++)
    if(params[i] != this->params[i].getCurrentEstimate())
      params[i] = this->params[i].getCurrentEstimate();

  for(map<int, ParamCalculatorShPtr>::const_iterator i = my_param_calcs.begin(); i != my_param_calcs.end(); ++i)
    params[i->first] = i->second->calculateParam(this);

  return 0;
}

//==============================================================
//  dumpConfiguration()
//==============================================================
void 
ParameterMgr::dumpConfiguration() const
{
  DebugCfg cfg(DebugCfg::COMPLETE);

  dump(cfg);
}

//=======================================================================
//  dump()
//=======================================================================
void 
ParameterMgr::dump(DebugCfg& debug) const
{
  for(vector<Parameter>::const_iterator p = params.begin(); p != params.end(); p++)
    p->dump(debug);
}

//=======================================================================
//  getGroupNames()
//=======================================================================
vector<string> 
ParameterMgr::getGroupNames() const
{
  vector<string> names(0);

  for(map<string, param_group>::const_iterator g = groups.begin(); g != groups.end(); g++)
    names.push_back(g->first);

  return names;
}

//=======================================================================
//  getParamBegin() NON CONST
//=======================================================================
ParameterIterator 
ParameterMgr::getParamBegin()
{
  return getParamBegin("ALL");
}

//=======================================================================
//  getParamBegin() CONST
//=======================================================================
ParameterConstIterator 
ParameterMgr::getParamBegin() const
{
  return getParamBegin("ALL");
}

//=======================================================================
//  getParamBegin() NON CONST
//=======================================================================
ParameterIterator 
ParameterMgr::getParamBegin(string group_name)
{
  return ParameterIterator(&params, getGroup(group_name).begin());
}

//=======================================================================
//  getParamBegin() CONST
//=======================================================================
ParameterConstIterator 
ParameterMgr::getParamBegin(string group_name) const
{
  return ParameterConstIterator(&params, getGroup(group_name).begin());
}

//=======================================================================
//  getParamEnd() NON CONST
//=======================================================================
ParameterIterator 
ParameterMgr::getParamEnd()
{
  return getParamEnd("ALL");
}

//=======================================================================
//  getParamEnd() CONST
//=======================================================================
ParameterConstIterator 
ParameterMgr::getParamEnd() const
{
  return getParamEnd("ALL");
}

//=======================================================================
//  getParamEnd() NON CONST
//=======================================================================
ParameterIterator 
ParameterMgr::getParamEnd(string group_name)
{
  return ParameterIterator(&params, getGroup(group_name).end());
}

//=======================================================================
//  getParamEnd() CONST
//=======================================================================
ParameterConstIterator 
ParameterMgr::getParamEnd(string group_name) const
{
  return ParameterConstIterator(&params, getGroup(group_name).end());
}

//=======================================================================
//  getParamCount()
//=======================================================================
int 
ParameterMgr::getParamCount() const
{
  return getParamCount("ALL");
}

//=======================================================================
//  getParamCount()
//=======================================================================
int 
ParameterMgr::getParamCount (string group_name) const
{
  return getGroup(group_name).size();
}

//=======================================================================
//  getEstimatedParamCount()
//=======================================================================
int 
ParameterMgr::getEstimatedParamCount() const
{
  int count = 0;

  for(ParameterConstIterator p = getParamBegin(); p != getParamEnd(); p++)
    count += p->isEstimated();

  return count;
}

//=======================================================================
//  getGroupCount()
//=======================================================================
int 
ParameterMgr::getGroupCount() const
{
  return groups.size();
}


NewSubmodelShPtr 
ParameterMgr::getSubmodel(const string& name)
{
  for(size_t i = 0; i < my_new_submodels.size(); ++i)
  {
    if(my_new_submodels[i]->getName() == name)
      return my_new_submodels[i];
  }

  SAGE_internal_error();

  // Note: never reached:
  return NewSubmodelShPtr();
}

NewSubmodelShCstPtr 
ParameterMgr::getSubmodel(const string& name) const
{
  for(size_t i = 0; i < my_new_submodels.size(); ++i)
  {
    if(my_new_submodels[i]->getName() == name)
      return my_new_submodels[i];
  }

  SAGE_internal_error();

  // Note: never reached:
  return NewSubmodelShCstPtr();
}

void
ParameterMgr::dumpParameterLookupTable() const
{
  cout << "\n\nParameter Lookup Table -" << endl;

  std::map<pair<string, string>, int>::const_iterator  p_iter     = param_lookup_table.begin();
  std::map<pair<string, string>, int>::const_iterator  p_end_iter = param_lookup_table.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    cout << "group " << p_iter->first.first  << "    "
         << "name "  << p_iter->first.second << "    "
         << "index " << p_iter->second       << endl << flush;
  }

  cout << "\n" << endl;  
}

// due to JA (July 09) for bypassing maximization
bool ParameterMgr::allfixed()
{
     bool all_fixed = true;
     vector<Parameter>::iterator itr;
   
      for (itr = params.begin(); itr != params.end(); itr++){ 
        if (itr->isEstimated()) {all_fixed = false; break;}
       }

      return all_fixed;
}


}
}
