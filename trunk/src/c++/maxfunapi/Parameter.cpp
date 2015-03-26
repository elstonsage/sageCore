//=======================================================================
//
//  File:	Parameter.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "maxfunapi/Parameter.h"

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//  ParamTypeEnum2str()
//=======================================================================
std::string
ParamTypeEnum2str(Parameter::ParamTypeEnum p)
{
  switch(p)
  {
    case Parameter::NO_PARAMTYPE                  : return "???";                           
    case Parameter::INDEPENDENT_FUNCTIONAL        : return "Independent functional";        
    case Parameter::INDEPENDENT                   : return "Independent";                   
    case Parameter::DEPENDENT                     : return "Dependent";                     
    case Parameter::FIXED                         : return "Fixed";                         
    case Parameter::IND_FUNC_FIXED_AT_BOUND       : return "Ind. func. fixed @ bnd";        
    case Parameter::IND_FIXED_AT_BOUND            : return "Ind. fixed @ bound";            
    case Parameter::IND_FUNC_FIXED_NEAR_BOUND     : return "Ind. func. fixed near bnd";     
    case Parameter::IND_FIXED_NEAR_BOUND          : return "Ind. fixed near bound";         
    case Parameter::IND_FUNC_FIXED_NOT_NEAR_BOUND : return "Ind. func. fixed not near bnd"; 
    case Parameter::IND_FIXED_NOT_NEAR_BOUND      : return "Ind. fixed not near bound";     
    default                                       : return "???";
  }
}

//=======================================================================
//  Parameter() CONSTRUCTOR #1
//=======================================================================
Parameter::Parameter()
{
  my_Name                      = "";
  my_NameAbbr                  = "";
  my_InitialEstimate           = 0.0;
  my_CurrentEstimate           = SAGE::QNAN;
  my_PreviousEstimate          = 0.0;
  my_FinalEstimate             = 0.0;
  my_GroupName                 = "";
  my_GroupIndex                = 0;
  my_InitialType               = NO_PARAMTYPE;
  my_IncludeInOutput           = true;
  my_IncludeDeriv              = true;
  my_IncludeStdError           = true;
  my_IncludepValue             = true;
  my_FinalType                 = NO_PARAMTYPE;
  my_LowerBound                = 0.0;
  my_UpperBound                = 0.0;
  my_InitialStepsizeFactor     = 0.1;
  my_StdErrorAvailable         = false;
  my_StdError                  = 0.0;
  my_DerivAvailable            = false;
  my_Deriv                     = 0.0;
  my_pValueAvailable           = false;
  my_pValue                    = 0.0;
  my_pvalue_opt                = TWO_SIDED;
  my_pvalue_mean               = 0.0;
  my_Index                     = 0;
  my_DerivIndex                = 0;
  my_VarAvailable              = false;
  my_VarIndex                  = 0;
  scoretest                    = false; // do not do score test
}

//=======================================================================
//  Parameter() COPY CONSTRUCTOR
//=======================================================================
Parameter::Parameter(const Parameter & other)
{
  copy(other);
}

//=======================================================================
//  ~Parameter() DESTRUCTOR
//=======================================================================
Parameter::~Parameter()
{}

//=======================================================================
//  Parameter() OPERATOR=
//=======================================================================
Parameter&
Parameter::operator=(const Parameter& other)
{
  if(&other != this)
  {
    copy(other);
  }

  return *this;
}

//=======================================================================
//  copy()
//=======================================================================
void
Parameter::copy(const Parameter& other)
{
  my_Name                      = other.my_Name;
  my_NameAbbr                  = other.my_NameAbbr;
  my_InitialEstimate           = other.my_InitialEstimate,
  my_PreviousEstimate          = other.my_PreviousEstimate;
  my_CurrentEstimate           = other.my_CurrentEstimate,
  my_FinalEstimate             = other.my_FinalEstimate;
  my_InitialType               = other.my_InitialType,
  my_IncludeInOutput           = other.my_IncludeInOutput;
  my_IncludeDeriv              = other.my_IncludeDeriv;
  my_IncludeStdError           = other.my_IncludeStdError;
  my_IncludepValue             = other.my_IncludepValue;
  my_FinalType                 = other.my_FinalType;
  my_GroupName                 = other.my_GroupName;
  my_GroupIndex                = other.my_GroupIndex;
  my_LowerBound                = other.my_LowerBound,
  my_UpperBound                = other.my_UpperBound;
  my_InitialStepsizeFactor     = other.my_InitialStepsizeFactor;
  my_StdErrorAvailable         = other.my_StdErrorAvailable;
  my_StdError                  = other.my_StdError;
  my_DerivAvailable            = other.my_DerivAvailable;
  my_Deriv                     = other.my_Deriv;
  my_pValueAvailable           = other.my_pValueAvailable;
  my_pValue                    = other.my_pValue;
  my_pvalue_opt                = other.my_pvalue_opt;
  my_pvalue_mean               = other.my_pvalue_mean;
  my_Index                     = other.my_Index;
  my_DerivIndex                = other.my_DerivIndex;
  my_VarAvailable              = other.my_VarAvailable;
  my_VarIndex                  = other.my_VarIndex;
  scoretest                    = other.scoretest;
}

//=======================================================================
//  isInBounds()
//=======================================================================
bool 
Parameter::isInBounds() const
{
  return ((my_CurrentEstimate >= my_LowerBound) && (my_CurrentEstimate <= my_UpperBound));
}

//=======================================================================
//  updateFinalOutput()
//=======================================================================
void
Parameter::updateFinalOutput(const Maxfun_Data & data, CovarianceMatrix& vcmatrix, int cov_matrix_status)
{
  // 0. Fetch final estimate, standard error, and final status:

	my_FinalEstimate =                data.param (my_Index);
	my_StdError      =                data.stde  (my_Index);
	my_FinalType     = (ParamTypeEnum)data.ist   (my_Index);

  // 1. Check that standard error is valid:

	if(finite(getStdError()) && getStdError() > 0)
	{
	  my_StdErrorAvailable = true;
	}
	else
	{
	  my_StdErrorAvailable = false;
	}

  // 2. Fetch derivative, if available:

	if(my_DerivIndex != (size_t)-1)
	{
	  my_Deriv          = data.g(my_DerivIndex);
	  my_DerivAvailable = true;
	}
	else
	{
	  my_DerivAvailable = false;
	}

  // 4. Calculate p-value:

        if(my_StdErrorAvailable == true)
        {
          calculatePValue();
        }
        else
        {
          my_pValueAvailable = false;
        }
}

//=======================================================================
//  calculatePValue()
//=======================================================================
void
Parameter::calculatePValue()
{
  if(my_pvalue_opt == ONE_SIDED)
    calculatePValueOneSided();
  else
    calculatePValueTwoSided();
}

//=======================================================================
//  calculatePValueOneSided()
//=======================================================================
void
Parameter::calculatePValueOneSided()
{
  if(my_StdErrorAvailable)
  {
    my_pValue = 1.0 - normal_cdf(abs((my_FinalEstimate - my_pvalue_mean) / my_StdError));

    if(!SAGE::isnan(my_pValue))
      my_pValueAvailable = true;
  }
}

//=======================================================================
//  calculatePValueTwoSided()
//=======================================================================
void
Parameter::calculatePValueTwoSided()
{
  if(my_StdErrorAvailable)
  {
    my_pValue = 2.0 * (1.0 - normal_cdf(abs((my_FinalEstimate - my_pvalue_mean) / my_StdError)));

    if(!SAGE::isnan(my_pValue))
      my_pValueAvailable = true;
  }
}

//=======================================================================
//  isFinalFixed()
//=======================================================================
bool
Parameter::isFinalFixed() const
{
  if(getFinalType() == INDEPENDENT || getFinalType() == INDEPENDENT_FUNCTIONAL || getFinalType() == DEPENDENT)
    return false;
  else
    return true;
}

//=======================================================================
//  dump()
//=======================================================================
void
Parameter::dump(DebugCfg& debug) const
{
  debug.getOutputStream() << left << setw(25) << my_GroupName << setw(20) << left << my_Name << " = "
           << setw(10) << left << std::fixed << my_CurrentEstimate
           << endl;
}

}
}
