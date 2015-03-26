#include "assoc/Model.h"

namespace SAGE  {
namespace ASSOC {

double 
ExpressionCalculator::operator() (size_t i, const SAMPLING::MemberDataSample & sample) const
{
  FUNC::Expression e;
    
  e.setExpression(my_expr);

  for(FUNC::Expression::VariableNames::const_iterator var = e.getVariableNames().begin(); var != e.getVariableNames().end(); ++var)
  {
    e.setVariable(*var, sample.getReplacedValue(i, "Covariates", *var));
  }
   
  return e.evaluate();
}

}
} 

