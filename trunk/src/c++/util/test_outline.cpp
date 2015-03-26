#include "util/OutlineCntr.h"
#include <iostream>

int main(int argc, char* argv[])
{
  SAGE::UTIL::OutlineCntr o;
  
  for(int i = 0; i < 3; ++i)
  {
    std::cout << (std::string)o << std::endl;
    
    o.increaseDepth();
    
    for(int j = 0; j < 3; ++j)
    {
      std::cout << (std::string)o << std::endl;

      o.increaseDepth();
      
      for(int k = 0; k < 3; ++k)
      {
        std::cout << (std::string)o << std::endl;

        o++;
      }

      o.decreaseDepth();

      o++;
    }
    
    o.decreaseDepth();

    o++;
  }  

  return 0;
}
