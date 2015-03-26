#ifndef ASSOC_MODEL_H
#define ASSOC_MODEL_H
//=======================================================================
///
///  File:	Model.h
///
///  Author:	Stephen Gross
///
///  Copyright 2002 R. C. Elston
//=======================================================================


#include <iostream>
#include <vector>
#include "LSF/parse_ops.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "fped/fped.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/maxfunapi.h"
#include "mfsubmodels/mfsubmodels.h"
#include "numerics/fmatrix.h"
#include "sampling/sampling.h"
#include "func/Expression.h"
#include <boost/cast.hpp>
#include "assoc/Datatypes.h"
#include "assoc/MemberCovariateCalculator.h"

namespace SAGE  {
namespace ASSOC {

class ResidTypeCalculator : public SAMPLING::TraitValueCalculator
{
  virtual double operator() (size_t i, const SAMPLING::MemberDataSample & sample) const
  {
    const FPED::Member & ind = sample.getIndividual(i);
                  
    bool is_founder   = !ind.parent1         (),
         has_children =  ind.offspring_count (),
         has_sibs     =  false;
                                              
    if(!is_founder)
    {
      if(ind.sibling_count())
      {
        has_sibs = true;
      }
      else if(ind.parent1()->mate_count() > 1 || ind.parent2()->mate_count() > 1)
      {
        has_sibs = (ind.parent1()->mate_count() == 1 ? ind.parent2() : ind.parent1())->offspring_count() > 1;
      }
    }
                                                                                                                                    
         if (has_children &&  is_founder)              return 0;
    else if (has_children && !is_founder && !has_sibs) return 1;
    else if (has_children && !is_founder &&  has_sibs) return 2;
    else if(!has_children && !is_founder && !has_sibs) return 3;
    else                                               return 4;
  }
};

class ExpressionCalculator : public SAMPLING::TraitValueCalculator
{
  public:
  
    typedef SAMPLING::TraitValueCalculator BaseType;
  
    explicit ExpressionCalculator(const std::string & expr) : BaseType(), my_expr(expr) { }
  
    ExpressionCalculator(const ExpressionCalculator & other) : BaseType(other), my_expr(other.my_expr) { }
  
    ExpressionCalculator & operator= (const ExpressionCalculator & other)
    {
      TraitValueCalculator::operator=(other);
    
      if(this != &other)
      {
        my_expr = other.my_expr;
      }
    
      return *this;
    }

    virtual double operator() (size_t i, const SAMPLING::MemberDataSample & sample) const;

  private:
  
    std::string my_expr;
};


}} /// End namespace

#endif
