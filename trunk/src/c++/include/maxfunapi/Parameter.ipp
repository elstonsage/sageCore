#ifndef PARAMETER_H
#include "maxfunapi/Parameter.h"
#endif

namespace SAGE   {
namespace MAXFUN {

//===================================================================
// Public accessors:
//===================================================================
    
inline pair<Parameter::PValueOptionEnum, double> Parameter::getPValueInf() const
{
  std::pair<Parameter::PValueOptionEnum, double> temp = make_pair(my_pvalue_opt, my_pvalue_mean);

  return temp;
}

inline bool Parameter::isEstimated() const { return (my_InitialType == INDEPENDENT) || (my_InitialType == INDEPENDENT_FUNCTIONAL); }

inline Parameter::ParamTypeEnum   Parameter::getInitialType  () const { return my_InitialType; }
inline Parameter::ParamTypeEnum   Parameter::getFinalType    () const { return my_FinalType;   }

inline const string      Parameter::getName                      () const { return my_Name;                      }
inline const string      Parameter::getNameAbbr                  () const { return my_NameAbbr;                  }
inline const string      Parameter::getGroupName                 () const { return my_GroupName;                 }

inline bool Parameter::estimateChanged() const { return getCurrentEstimate() != getPreviousEstimate(); }

//inline bool        Parameter::isInUse                      () const { return my_InUse;                     }
inline double      Parameter::getInitialEstimate           () const { return my_InitialEstimate;           }
inline double      Parameter::getPreviousEstimate          () const { return my_PreviousEstimate;          }
inline double      Parameter::getCurrentEstimate           () const { return my_CurrentEstimate;           }
inline double      Parameter::operator()                   () const { return my_CurrentEstimate;           }
inline double      Parameter::getFinalEstimate             () const { return my_FinalEstimate;             }
inline size_t      Parameter::getGroupIndex                () const { return my_GroupIndex;                }
inline double      Parameter::getLowerBound                () const { return my_LowerBound;                }
inline double      Parameter::getUpperBound                () const { return my_UpperBound;                }
inline double      Parameter::getInitialStepsizeFactor     () const { return my_InitialStepsizeFactor;     }
inline bool        Parameter::isStdErrorAvailable          () const { return my_StdErrorAvailable;         }
inline double      Parameter::getStdError                  () const { return my_StdError;                  }
inline bool        Parameter::isDerivAvailable             () const { return my_DerivAvailable;            }
inline double      Parameter::getDeriv                     () const { return my_Deriv;                     }
inline bool        Parameter::isPValueAvailable            () const { return my_pValueAvailable;           }
inline double      Parameter::getPValue                    () const { return my_pValue;                    }
inline size_t      Parameter::getIndex                     () const { return my_Index;                     }
inline size_t      Parameter::getDerivIndex                () const { return my_DerivIndex;                }
inline bool        Parameter::isVarAvailable               () const { return my_VarAvailable;              }
inline size_t      Parameter::getVarIndex                  () const { return my_VarIndex;                  }
inline bool        Parameter::getIncludeInOutput           () const { return my_IncludeInOutput;           }
inline bool        Parameter::getIncludeDeriv              () const { return my_IncludeDeriv;              }
inline bool        Parameter::getIncludeStdError           () const { return my_IncludeStdError;           }
inline bool        Parameter::getIncludePValue             () const { return my_IncludepValue;             }

//===================================================================
// Public mutators:
//===================================================================

inline int Parameter::setPValueInf(PValueOptionEnum opt, double mean)
{
  my_pvalue_opt  = opt;
  my_pvalue_mean = mean;

  return 0;
}
                
inline int Parameter::setFinalType (Parameter::ParamTypeEnum x) { my_FinalType = x; return 0; }
// due to JA for facilitating the score test

inline int Parameter::setInitialType (Parameter::ParamTypeEnum x) { my_InitialType = x; return 0; }

inline double & Parameter::operator() () { return my_CurrentEstimate; }

inline int Parameter::setName            (string x) { my_Name             = x; return 0; }
inline int Parameter::setNameAbbr        (string x) { my_NameAbbr         = x; return 0; }
//inline int Parameter::setInUse           (bool   x) { my_InUse            = x; return 0; }
inline int Parameter::setCurrentEstimate (double x) { my_PreviousEstimate = my_CurrentEstimate;
                                                      my_CurrentEstimate  = x; return 0; }
inline int Parameter::setFinalEstimate   (double x) { my_FinalEstimate    = x; return 0; }
inline int Parameter::setStdError        (double x) { my_StdError         = x; return 0; }
inline int Parameter::setDeriv           (double x) { my_Deriv            = x; return 0; }
inline int Parameter::setGroupName       (string x) { my_GroupName        = x; return 0; }
inline int Parameter::setGroupIndex      (size_t x) { my_GroupIndex       = x; return 0; }
inline int Parameter::setLowerBound      (double x) { my_LowerBound       = x; return 0; }
inline int Parameter::setUpperBound      (double x) { my_UpperBound       = x; return 0; }
inline int Parameter::setInitialStepsizeFactor(double x)
                                                    { my_InitialStepsizeFactor = x; return 0;}
inline int Parameter::setIndex           (size_t x) { my_Index            = x; return 0; }
inline int Parameter::setDerivIndex      (size_t x) { my_DerivIndex       = x; return 0; }
inline int Parameter::setVarIndex        (size_t x) { my_VarIndex         = x; return 0; }
inline int Parameter::setIncludeInOutput (bool   x) { my_IncludeInOutput  = x; return 0; }
inline int Parameter::setIncludeDeriv    (bool   x) { my_IncludeDeriv     = x; return 0; }
inline int Parameter::setIncludeStdError (bool   x) { my_IncludeStdError  = x; return 0; }
inline int Parameter::setIncludePValue   (bool   x) { my_IncludepValue    = x; return 0; }

inline int
Parameter::setInitialEstimate(double initial_estimate)
{ 
  my_InitialEstimate = initial_estimate;

  return 0;
}

//===============================================================================================================
//===============================================================================================================
//===============================================================================================================
//
//			NON-CONST ITERATOR
//
//===============================================================================================================
//===============================================================================================================
//===============================================================================================================

//===========================================================================
// CONSTRUCTOR - Protected
//===========================================================================
inline 
ParameterIterator::ParameterIterator(
	param_vector          * _masterlist, 
	param_group::iterator   _iter) 
	:
	masterlist (_masterlist),
	iter       (_iter)
{}

//===========================================================================
// COPY CONSTRUCTOR
//===========================================================================
inline
ParameterIterator::ParameterIterator(const ParameterIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// operator* ()
//===========================================================================
inline ParameterIterator::reference
ParameterIterator::operator*()
{ 
  return (*masterlist)[*iter];
}

//===========================================================================
// operator[] (...)
//===========================================================================
inline ParameterIterator::reference
ParameterIterator::operator[](difference_type i)
{ 
  return *(*this + i); 
}
        
//===========================================================================
// operator-> ()
//===========================================================================
inline ParameterIterator::pointer 
ParameterIterator::operator->()
{
  return &(*masterlist)[*iter];
}

//===========================================================================
// operator++ () POST INCREMENT
//===========================================================================
inline ParameterIterator::iterator
ParameterIterator::operator++()
{ 
  increment();
  return *this;
}
  
//===========================================================================
// operator++ (...) PRE INCREMENT
//===========================================================================
inline ParameterIterator::iterator
ParameterIterator::operator++(int)
{ 
  iterator tmp = *this; 
  increment(); 
  return tmp; 
}

//===========================================================================
// operator-- () POST DECREMENT
//===========================================================================
inline ParameterIterator::iterator
ParameterIterator::operator--()
{ 
  decrement();
  return *this;
}

//===========================================================================
// operator-- (...) PRE DECREMENT
//===========================================================================
inline ParameterIterator::iterator
ParameterIterator::operator--(int)
{ 
  iterator tmp = *this; 
  decrement(); 
  return tmp; 
}

//===========================================================================
// operator- (...)
//===========================================================================
inline ParameterIterator::difference_type 
ParameterIterator::operator-(const iterator & other) const
{
  return iter - other.iter;
}

//===========================================================================
// operator== (...)
//===========================================================================
inline bool
ParameterIterator::operator== (const iterator & other) const
{
  return iter == other.iter;
}

//===========================================================================
// operator!= (...)
//===========================================================================
inline bool 
ParameterIterator::operator!= (const iterator & other) const
{
  return !(*this == other);
}

//===========================================================================
// operator< (...)
//===========================================================================
inline bool 
ParameterIterator::operator< (const iterator & other) const
{
  return iter < other.iter;
}

//===========================================================================
// increment()
//===========================================================================
inline void 
ParameterIterator::increment()
{
  iter++;
}

//===========================================================================
// decrement()
//===========================================================================
inline void 
ParameterIterator::decrement()
{
  iter--;
}

//===============================================================================================================
//===============================================================================================================
//===============================================================================================================
//
//			CONST ITERATOR
//
//===============================================================================================================
//===============================================================================================================
//===============================================================================================================

//===========================================================================
// CONSTRUCTOR - Protected
//===========================================================================
inline 
ParameterConstIterator::ParameterConstIterator(
	const param_vector          * _masterlist, 
	param_group::const_iterator   _iter) 
	:
	masterlist (_masterlist),
	iter       (_iter)
{}

//===========================================================================
// COPY CONSTRUCTOR #1
//===========================================================================
inline
ParameterConstIterator::ParameterConstIterator(const ParameterConstIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// COPY CONSTRUCTOR #2
//===========================================================================
inline
ParameterConstIterator::ParameterConstIterator(const ParameterIterator & other) :
	masterlist (other.masterlist),
	iter       (other.iter)
{}

//===========================================================================
// operator* ()
//===========================================================================
inline ParameterConstIterator::const_reference
ParameterConstIterator::operator*() const
{ 
  return (*masterlist)[*iter];
}

//===========================================================================
// operator[] (...)
//===========================================================================
inline ParameterConstIterator::const_reference
ParameterConstIterator::operator[](difference_type i) const
{ 
  return *(*this + i); 
}
        
//===========================================================================
// operator-> ()
//===========================================================================
inline ParameterConstIterator::const_pointer 
ParameterConstIterator::operator->()
{
  return &(*masterlist)[*iter];
}

//===========================================================================
// operator++ () POST INCREMENT
//===========================================================================
inline ParameterConstIterator::const_iterator
ParameterConstIterator::operator++()
{ 
  increment();
  return *this;
}
  
//===========================================================================
// operator++ (...) PRE INCREMENT
//===========================================================================
inline ParameterConstIterator::const_iterator
ParameterConstIterator::operator++(int)
{ 
  const_iterator tmp = *this; 
  increment(); 
  return tmp; 
}

//===========================================================================
// operator-- () POST DECREMENT
//===========================================================================
inline ParameterConstIterator::const_iterator
ParameterConstIterator::operator--()
{ 
  decrement();
  return *this;
}

//===========================================================================
// operator-- (...) PRE DECREMENT
//===========================================================================
inline ParameterConstIterator::const_iterator
ParameterConstIterator::operator--(int)
{ 
  const_iterator tmp = *this; 
  decrement(); 
  return tmp; 
}

//===========================================================================
// operator- (...)
//===========================================================================
inline ParameterConstIterator::difference_type 
ParameterConstIterator::operator-(const const_iterator & other) const
{
  return iter - other.iter;
}

//===========================================================================
// operator== (...)
//===========================================================================
inline bool
ParameterConstIterator::operator== (const const_iterator & other) const
{
  return iter == other.iter;
}

//===========================================================================
// operator!= (...)
//===========================================================================
inline bool 
ParameterConstIterator::operator!= (const const_iterator & other) const
{
  return !(*this == other);
}

//===========================================================================
// operator< (...)
//===========================================================================
inline bool 
ParameterConstIterator::operator< (const const_iterator & other) const
{
  return iter < other.iter;
}

//===========================================================================
// increment()
//===========================================================================
inline void 
ParameterConstIterator::increment()
{
  iter++;
}

//===========================================================================
// decrement()
//===========================================================================
inline void 
ParameterConstIterator::decrement()
{
  iter--;
}

}} // End namespace
