#include "maxfunapi/Submodel.h"
#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//  MaxfunSubModel() CONSTRUCTOR
//=======================================================================
Submodel::Submodel(cerrorstream& errors) : my_errors (errors), my_info   (NULL)
{ 
  my_parameters.clear();
}

//=======================================================================
//  MaxfunSubModel() COPY CONSTRUCTOR
//=======================================================================
Submodel::Submodel(const Submodel& other) 
      : my_parameters(other.my_parameters), my_errors(other.my_errors), my_info(NULL)
{}

//=======================================================================
//  operator=()
//=======================================================================
Submodel& 
Submodel::operator=(const Submodel& other)
{
  if(this != &other)
  {
    my_info       = NULL;
    my_errors     = other.my_errors;
    my_parameters = other.my_parameters;
  }

  return *this;
}  

//=======================================================================
//  isLinked()
//=======================================================================
bool
Submodel::isLinked() const
{
  return (my_info != NULL);
}

//=======================================================================
//  linkToParameterMgr()
//=======================================================================
void 
Submodel::linkToParameterMgr(ParameterMgr* mi)
{
  my_info = mi;
}

//============================================================================
// addParametersToParameterMgr()
//============================================================================
int
Submodel::addParametersToParameterMgr()
{
  for(size_t i = 0; i < my_parameters.size(); ++i)
    my_info->addParameter(my_parameters[i]);

  return 0;
}

//============================================================================
// getParam() CONST and NON-CONST
//
//============================================================================
double & Submodel::getParam(int local_id)       { return my_info->operator()(my_parameters[local_id].index); } 
double   Submodel::getParam(int local_id) const { return my_info->operator()(my_parameters[local_id].index); } 

//=======================================================================
//  unlinkFromParameterMgr()
//=======================================================================
void 
Submodel::unlinkFromParameterMgr()
{
  my_info = NULL;
}

//=========================================================
//  finalizeConfiguration()
//=========================================================
int 
Submodel::finalizeConfiguration()
{
  return 0;
}

//=========================================================
//  setAdvancedParameterOptions()
//=========================================================
void 
Submodel::setAdvancedParameterOptions()
{}

//=========================================================
//  ~Submodel() DESTRUCTOR
//=========================================================
Submodel::~Submodel()
{
  if(my_info)
  {
    ParameterMgr* tmp = my_info;

    my_info = NULL;

    tmp->removeSubmodel(this);
  }
}

} // End namespace MAXFUN
} // End namespace SAGE
