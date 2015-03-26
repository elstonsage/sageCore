#include <iostream>
#include <fstream>
#include <string>
#include "numerics/sinfo.h"

void print_info(const SAGE::SampleInfo& info)
{
  std::cout << "count    = " << info.count()            << "     \tsum    = " << info.sum()  << std::endl;
  std::cout << "mean     = " << info.mean()                  << "\tsum^2  = " << info.sum2() << std::endl;
  std::cout << "variance = " << info.variance()              << "\tsum^3  = " << info.sum3() << std::endl;
  std::cout << "std.dev  = " << info.standard_deviation()    << "\tsum^4  = " << info.sum4() << std::endl;
  std::cout << "std.err  = " << info.standard_error()        << std::endl;
  std::cout << "ceff.var = " << info.variation_coefficient() << std::endl;
  std::cout << "skewness = " << info.skewness() << std::endl;
  std::cout << "kutrosis = " << info.kurtosis() << std::endl;
  std::cout << "range    = " << info.min() << " .. " << info.max() << std::endl;
}

void process(std::istream& in)
{
  if( in.bad() )
  {
    std::cerr << "Error: Cannot read input file" << std::endl;
    return;
  }

  double d;
  SAGE::SampleInfo info;
  info.sample_adjustment(1);

  while(!in.eof())
  {

    in >> d;

    if(!in || !finite(d) )
    {
      if( !in.eof() && in.fail() )
        in.clear();
      continue;
    }     

    info += d;
  }

  print_info(info);
}

int main(int argc, char* argv[])
{
  if(argc < 2)
    process(std::cin);
  else
  {
    std::ifstream infile( argv[1] );
    process(infile);
  }
  return 0;
}
