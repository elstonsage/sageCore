#include "util/StringUtils.h"
#include <iostream>

int main(int argc, char* argv[])
{
  std::cout << "Testing multi delimit split:" << std::endl;

  std::string delimiter_list = " ,;";

  std::cout << "Delimiter list = <" << delimiter_list << ">" << std::endl;

  std::vector<std::string> input_lines;

  input_lines.push_back("hello");
  input_lines.push_back("hello there this is a big picture");
  input_lines.push_back("hello,there,this,;is,a,big,picture");
  input_lines.push_back("hello,there this;  ,is a,big ;;picture");
  input_lines.push_back("hello, there , , ,;;;;;; this is a big ,,,,,,    picture");

  std::vector<std::string> results;

  for(size_t i = 0; i < input_lines.size(); ++i)
  {
    std::cout << input_lines[i] << " --> token count = " 
              << SAGE::UTIL::StringUtils::splitMultiDelimitedString(input_lines[i], delimiter_list, results);

    std::cout << " (";

    for(size_t j = 0; j < results.size(); ++j)
      std::cout << (j?",":"") << results[j];

    std::cout << ")" << std::endl;
  }

  std::cout << SAGE::UTIL::StringUtils::lineWrapString("This is a very special string that needs to be split", 10) << std::endl;

  SAGE::UTIL::StringUtils::splitMultiDelimitedString(" , ,  token,token", " ,", results);
  
  for(size_t i = 0; i < results.size(); ++i)
    std::cout << results[i] << ",";
    
  std::cout << std::endl;

  return 0;
}
