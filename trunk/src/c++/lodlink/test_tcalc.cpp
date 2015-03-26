//============================================================================
// File:      test_tcalc.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 10/16/2                                                   
//                                                                          
// Notes:     Tests the lodlink trans_calculator class.
//    
//                                                                      
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <fstream>
#include "lodlink/trans_calculator.h"

using namespace std;
using namespace SAGE;
using namespace LODLINK;

void  check_transition(ofstream& out);

int main(int argc, char* argv[])
{

  ofstream out;
  out.open("testtcalc.out");
  
  if(! out)
  {
    cout << "Cannot open output file: testtcalc.out.  Exiting ..." << endl;
    exit(EXIT_FAILURE);
  }

  check_transition(out);

  exit(EXIT_SUCCESS);
}

void  print_header(ofstream& out, int test);
void  print_trio(ofstream& out, const joint_genotype& mom,
                                const joint_genotype& dad,
                                const joint_genotype& kid );
void  print_result(ofstream& out, int test,
                   const joint_genotype& mom,
                   const joint_genotype& dad,
                   const joint_genotype& kid,
                   const mle_sub_model& mle,
                   double prob               );

void
check_transition(ofstream& out)
{
  MLOCUS::genotype_model  tgm;
  MLOCUS::genotype_model  mgm;

  // Trait alleles.
  tgm.add_allele("a", .3);
  tgm.add_allele("b", .3);
  tgm.add_allele("c", .4);
  
  MLOCUS::allele  a = tgm.get_allele("a");
  MLOCUS::allele  b = tgm.get_allele("b");
  MLOCUS::allele  c = tgm.get_allele("c");
  
  MLOCUS::phased_genotype   aa = tgm.get_phased_genotype(a, a);
  MLOCUS::phased_genotype   bb = tgm.get_phased_genotype(b, b);
  MLOCUS::phased_genotype   cc = tgm.get_phased_genotype(c, c);
  MLOCUS::phased_genotype   ab = tgm.get_phased_genotype(a, b);
  MLOCUS::phased_genotype   ac = tgm.get_phased_genotype(a, c);
  MLOCUS::phased_genotype   ba = tgm.get_phased_genotype(b, a);
  MLOCUS::phased_genotype   bc = tgm.get_phased_genotype(b, c);
  MLOCUS::phased_genotype   ca = tgm.get_phased_genotype(c, a);
  MLOCUS::phased_genotype   cb = tgm.get_phased_genotype(c, b);
  
  // Marker alleles.
  mgm.add_allele("x", .3);
  mgm.add_allele("y", .3);
  mgm.add_allele("z", .4);
  
  MLOCUS::allele  x = mgm.get_allele("x");
  MLOCUS::allele  y = mgm.get_allele("y");
  MLOCUS::allele  z = mgm.get_allele("z");
  
  MLOCUS::phased_genotype   xx = mgm.get_phased_genotype(x, x);
  MLOCUS::phased_genotype   yy = mgm.get_phased_genotype(y, y);
  MLOCUS::phased_genotype   zz = mgm.get_phased_genotype(z, z);
  MLOCUS::phased_genotype   xy = mgm.get_phased_genotype(x, y);
  MLOCUS::phased_genotype   xz = mgm.get_phased_genotype(x, z);
  MLOCUS::phased_genotype   yx = mgm.get_phased_genotype(y, x);
  MLOCUS::phased_genotype   yz = mgm.get_phased_genotype(y, z);
  MLOCUS::phased_genotype   zx = mgm.get_phased_genotype(z, x);
  MLOCUS::phased_genotype   zy = mgm.get_phased_genotype(z, y);
  
  // Test 1.
  joint_genotype    mom(aa, xx);
  joint_genotype    dad(ba, yx);
  joint_genotype    kid(ab, xx);
  mle_sub_model     mle;
  trans_calculator  inst(mle);
  double            prob = inst.transition(mom, dad, kid);
  print_result(out, 1, mom, dad, kid, mle, prob);

  // Test 2.
  mom.tg = ba;
  mom.mg = xy;
  dad.tg = ba;
  dad.mg = xy;
  kid.tg = aa;
  kid.mg = xx;
  mle.set(true, false);
  mle.set_male_theta(.2);
  mle.set_female_theta(.3);
  prob = inst.transition(mom, dad, kid);
  print_result(out, 2, mom, dad, kid, mle, prob);  
  
  // Test 3.
  kid.mg = xz;
  prob = inst.transition(mom, dad, kid);
  print_result(out, 3, mom, dad, kid, mle, prob);
  
  // Test 4.
  mom.tg = ba;
  mom.mg = xy;
  dad.tg = ba;
  dad.mg = yx;
  kid.tg = ab;
  kid.mg = xy;
  prob = inst.transition(mom, dad, kid);
  print_result(out, 4, mom, dad, kid, mle, prob);  
}

void
print_header(ofstream& out, int test)
{
  out << "\n----------   test " << test << "   ----------\n" << endl;
}

void  
print_result(ofstream& out, int test,
             const joint_genotype& mom,
             const joint_genotype& dad,
             const joint_genotype& kid,
             const mle_sub_model& mle,
             double prob               )
{
  print_header(out, test);
  print_trio(out, mom, dad, kid);
  out << mle;
  out << "transition probability:  " << prob << endl;
}

void
print_trio(ofstream& out, const joint_genotype& mom,
                          const joint_genotype& dad,
                          const joint_genotype& kid )
{
  out << "mom: " << mom << "\n"
      << "dad: " << dad << "\n"
      << "kid: " << kid << "\n"
      << endl; 
}

