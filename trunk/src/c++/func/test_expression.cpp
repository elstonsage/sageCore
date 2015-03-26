#include "func/Expression.h"

int main()
{
  SAGE::FUNC::Expression e;
  
  e.dump(); 
  
  e.setExpression("x*y*sqrt(x)");

  e.setVariable("x", 5.0);
  e.setVariable("y", 5.0);
 
  e.dump();

  std::cout << e.evaluate() << std::endl;

  return 0;
}
