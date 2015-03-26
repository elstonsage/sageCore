#include "util/RegexUtils.h"
#include <iostream>

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Usage: test_regex pattern input_line" << std::endl;
    exit(0);
  }

  std::vector<std::string> f = SAGE::UTIL::RegexUtils::matchPattern(argv[1], argv[2]);

  for(size_t i = 0; i < f.size(); ++i)
    std::cout << f[i] << std::endl;

  return 0;
}
