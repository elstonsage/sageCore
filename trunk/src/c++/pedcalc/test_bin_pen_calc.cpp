#include "pedcalc/binary_penetrance_calculator.h"
#include "output/Output.h"
#include <iostream>

using namespace SAGE;
using namespace PED_CALC;
using namespace OUTPUT;
using namespace std;

// Realistically, this isn't the kind of parameters we'll use, but we're testing
// the numerics, so it doesn't matter.
double suscept(double member, double genotype)
{
  return member * genotype;
}

bool affection1(double member)
{
  return member < 0.5;
}
bool affection2(double member)
{
  return member > 0.5;
}

OUTPUT::Table generate_test_output
  (const BinaryPenetranceCalculator<double, double>& bin_calc,
   bool (*affection)(double),
   const std::string& name)
{
  OUTPUT::Table output(name);
  
  output << TableColumn("Member")
         << TableColumn("Genotype")
         << TableColumn("Suscept.")
         << TableColumn("Aff.")
         << TableColumn("Penetrance");
  
  // "Member" iteration
  for(double i = 0; i < 1.0; i += 0.1)
  {
    // "Genotype" iteration
    for(double j = -2; j < 2.0; j+= 0.2)
    {
      TableRow tr;
      
      tr << i << j << suscept(i,j) << affection(i) << bin_calc(i,j);
      
      output << tr;
    }
  }

  return output;
}
int main()
{
  BinaryPenetranceCalculator<double, double> bin_calc1(&suscept, &affection1);
  BinaryPenetranceCalculator<double, double> bin_calc2(&suscept, &affection2);

  BinaryPenetranceCalculator<double, double> bin_calc_copy(bin_calc1);
  
  cout << generate_test_output(bin_calc1,     affection1, "Affection1")
       << generate_test_output(bin_calc2,     affection2, "Affection2")
       << generate_test_output(bin_calc_copy, affection1, "Copy Constructor");
       
  bin_calc_copy = bin_calc2;
 
  cout << generate_test_output(bin_calc_copy,     affection2, "Copy Operator");
}
