#include "containers/AnyVector.h"

int main()
{
  SAGE::AnyVector<int, double, char> t;
  
  t.push_back(9);
  t.push_back(9.9);
  t.push_back('9');
  t.push_back(8);
  t.push_back(8.9);
  t.push_back('8');
  t.push_back(7);
  t.push_back(7.9);
  t.push_back('7');

  std::cout << "Ints: ";

  for(SAGE::AnyVectorItrs::ConstIterator<int> i = t.begin<int>(); i != t.end<int>(); ++i)
    std::cout << *i << " ";
    
  std::cout << std::endl;

  std::cout << "Doubles: ";

  for(SAGE::AnyVectorItrs::ConstIterator<double> i = t.begin<double>(); i != t.end<double>(); ++i)
    std::cout << *i << " ";
    
  std::cout << std::endl;

  std::cout << "Chars: ";

  for(SAGE::AnyVectorItrs::ConstIterator<char> i = t.begin<char>(); i != t.end<char>(); ++i)
    std::cout << *i << " ";
    
  std::cout << std::endl;

  for(size_t i = 0; i < t.size(); ++i)
  {
    if(t.isType<int>(i))
    {
      std::cout << "int " << t.getAbs<int>(i) << std::endl;
    }
    else if(t.isType<double>(i))
    {
      std::cout << "double " << t.getAbs<double>(i) << std::endl;
    }
    else if(t.isType<char>(i))
    {
      std::cout << "char " << t.getAbs<char>(i) << std::endl;
    }
  }

  return 0;
}
