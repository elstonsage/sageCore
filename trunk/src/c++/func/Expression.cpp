#include "func/Expression.h"
#include "func/PythonInterface.h"
#include "globals/SAGEConstants.h"
#include "error/internal_error.h"

namespace SAGE {
namespace FUNC {

Expression::Expression()
{
  my_expression      = "";
  my_variable_names  . clear();
  my_variable_values . clear();
  
  setExpression("0");
}

Expression::Expression(const Expression & other) : 
  my_expression      (other.my_expression),
  my_variable_names  (other.my_variable_names),
  my_variable_values (other.my_variable_values)
{ }

Expression& 
Expression::operator= (const Expression & other)
{
  if(this != &other)
  {
    my_expression      = other.my_expression;
    my_variable_names  = other.my_variable_names;
    my_variable_values = other.my_variable_values;
  }
  
  return *this;
}
            

bool 
Expression::setExpression(const std::string & expression)
{
  PythonInterface calculator(sage_cerr);

  if(calculator.compile(expression))
  {
    try
    {
      my_expression     = expression;
      my_variable_names . clear();

      VariableNames names = calculator.getNameList(my_expression);
      
      for(size_t i = 0; i < names.size(); ++i)
        if(!calculator.isNameInEnvironment(names[i]))
          my_variable_names.push_back(names[i]);

      my_variable_values.resize(my_variable_names.size(), SAGE::QNAN);

      return true;
    }
    catch(const PythonInterface::PythonException & e)
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

const Expression::VariableNames &
Expression::getVariableNames() const
{
  return my_variable_names;
}
    
bool 
Expression::setVariable(const std::string & name, double value)
{
  for(size_t i = 0; i < my_variable_names.size(); ++i)
  {
    if(my_variable_names[i] == name)
    {
      my_variable_values[i] = value;
      return true;
    }
  }
  
  return false;
}

double 
Expression::evaluate() const
{
  /// Create the interface:
  PythonInterface calculator(sage_cerr);
  
  // Compile the expression:
  calculator.compile(my_expression);
    
  // Set the variables:
  for(size_t i = 0; i < my_variable_names.size(); ++i)
    calculator.addDoubleToEnvironment(my_variable_names[i], my_variable_values[i]);
    
  // Calculate it!
  return calculator.run(my_expression, 30);
}

void
Expression::dump() const
{
  std::cout << "Expression:" << std::endl
            << "  my_expression = " << my_expression << std::endl
            << "  variables = ";
            
  for(size_t i = 0; i < my_variable_names.size(); ++i)
  {
    std::cout << my_variable_names[i] << " (";
    
    if(SAGE::isnan(my_variable_values[i]))
      std::cout << "unassigned";
    else
      std::cout << my_variable_values[i];
      
    std::cout << ") ";
  }
  
  std::cout << std::endl;
}
    
} // End namespace FUNC
} // End namespace SAGE
