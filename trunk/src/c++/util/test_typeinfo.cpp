#include "util/TypeInfo.h"
#include <iostream>

int main(int argc, char* argv[])
{
  typedef SAGE::UTIL::TypeInfo T;
  
  std::cout << "int == float? " << (T::create<int>() == T::create<float>()) << std::endl;
  std::cout << "int != float? " << (T::create<int>() != T::create<float>()) << std::endl;
  std::cout << "int <  float? " << (T::create<int>() <  T::create<float>()) << std::endl;
  std::cout << "int <= float? " << (T::create<int>() <= T::create<float>()) << std::endl;
  std::cout << "int >  float? " << (T::create<int>() >  T::create<float>()) << std::endl;
  std::cout << "int >= float? " << (T::create<int>() >= T::create<float>()) << std::endl;
  
  return 0;
}
