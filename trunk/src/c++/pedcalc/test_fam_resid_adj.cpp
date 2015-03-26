#include "pedcalc/fra_test_util.h"
#include "pedcalc/fam_resid_adj.h"
#include "mped/mp.h"
#include "boost/iterator/counting_iterator.hpp"
#include <iostream>
#include <sstream>

using namespace SAGE;
using namespace MPED;

enum Genotype { gAA, gAB, gBB, gEnd };

namespace test1
{

double residual(PED_CALC::FraResidType)
{
  return 0.1;
}
bool affection(const member_base& m)
{
  return m.is_male();
}
double susceptibility(const member_base&, Genotype g)
{
  return (g == gAA) ? -5 : 5;
}

double transmission(const member_base&, Genotype mother_geno,
                                        Genotype father_geno,
                                        Genotype child_geno)
{
  int mBcount =     (int) mother_geno;
  int mAcount = 2 - (int) mother_geno;

  int fBcount =     (int) father_geno;
  int fAcount = 2 - (int) father_geno;
  
  if      (child_geno == gAA) return 0.25 * mAcount * fAcount;
  else if (child_geno == gAB) return 0.25 * (mAcount * fBcount +  mBcount * fAcount);
  else                        return 0.25 * mBcount * fBcount;
}

void run()
{
  OUTPUT::Section test("Test 1 : 3 person test");

  pedigree_base p("1");
  
  p.add_member("1", SEX_MALE);
  p.add_member("2", SEX_FEMALE);
  
  p.add_member("3", SEX_FEMALE);

  p.add_lineage("3", "1", "2");

  p.build();
  p.freeze();
  
  PED_CALC::ExactFamResidAdj<Genotype, multipedigree_base>
      exact_adj
        (&residual,
         &affection,
         PED_CALC::BinaryPenetranceCalculator<member_base, Genotype>
             (&susceptibility, &affection));

  test << PED_CALC::generate_fra_test_output
              (exact_adj, p, boost::counting_iterator<int>(gAA),
                             boost::counting_iterator<int>(gEnd));
                             
 
                       
  PED_CALC::ApproximateFamResidAdj<Genotype, multipedigree_base, 
                                   boost::counting_iterator<int> >
      approx_adj
        (&residual,
         &affection,
         &susceptibility,
         &transmission,
         boost::counting_iterator<int>(gAA),
         boost::counting_iterator<int>(gEnd));
         
  test << PED_CALC::generate_fra_test_output
              (approx_adj, p, boost::counting_iterator<int>(gAA),
                              boost::counting_iterator<int>(gEnd));
 
  std::cout << test;
}

}

namespace test2
{

double residual(PED_CALC::FraResidType)
{
  return 0.1;
}
bool affection(const member_base& m)
{
  return m.is_male();
}
double susceptibility(const member_base&, Genotype g)
{
  return (g == gAA) ? -5 : 5;
}

double transmission(const member_base&, Genotype mother_geno,
                                        Genotype father_geno,
                                        Genotype child_geno)
{
  int mBcount =     (int) mother_geno;
  int mAcount = 2 - (int) mother_geno;

  int fBcount =     (int) father_geno;
  int fAcount = 2 - (int) father_geno;
  
  if      (child_geno == gAA) return 0.25 * mAcount * fAcount;
  else if (child_geno == gAB) return 0.25 * (mAcount * fBcount +  mBcount * fAcount);
  else                        return 0.25 * mBcount * fBcount;
}

void run()
{
  OUTPUT::Section test("Test 2 : 5 person pedigree (3 sibs) test");

  pedigree_base p("1");
  
  p.add_member("1", SEX_MALE);
  p.add_member("2", SEX_FEMALE);
  
  p.add_member("3", SEX_FEMALE);
  p.add_member("4", SEX_MALE);
  p.add_member("5", SEX_FEMALE);

  p.add_lineage("3", "1", "2");
  p.add_lineage("4", "1", "2");
  p.add_lineage("5", "1", "2");

  p.build();
  p.freeze();
  
  PED_CALC::ExactFamResidAdj<Genotype, multipedigree_base>
      exact_adj
        (&residual,
         &affection,
         PED_CALC::BinaryPenetranceCalculator<member_base, Genotype>
             (&susceptibility, &affection));

  test << PED_CALC::generate_fra_test_output
              (exact_adj, p, boost::counting_iterator<int>(gAA),
                             boost::counting_iterator<int>(gEnd));
                             
 
                       
  PED_CALC::ApproximateFamResidAdj<Genotype, multipedigree_base, 
                                   boost::counting_iterator<int> >
      approx_adj
        (&residual,
         &affection,
         &susceptibility,
         &transmission,
         boost::counting_iterator<int>(gAA),
         boost::counting_iterator<int>(gEnd));
         
  test << PED_CALC::generate_fra_test_output
              (approx_adj, p, boost::counting_iterator<int>(gAA),
                              boost::counting_iterator<int>(gEnd));
 
  std::cout << test;
}

}

int main()
{
  test1::run();
  test2::run();
  
  return 0;
}
