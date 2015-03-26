//===========================================================================
//
//  File:	polygenic_transition_calculator.cpp
//
//  Author:	Stephen Gross
//
//  Purpose:	Calculate polygenic transition values.
//
//  Copyright (c) 2001, R. C. Elston
//
//===========================================================================

#include "segreg/polygenic_transition_calculator.h"

using namespace std;
namespace SAGE {
namespace SEGREG {


//============================================================================
// Constructor
//============================================================================
polygenic_transition_calculator::polygenic_transition_calculator(const model & mod_param) :
	mod(mod_param)
{
  // Construct probs vector:

  vector<double>                   v1;
  vector<vector<double> >          v2;
  vector<vector<vector<double> > > v3;

  double QNAN = std::numeric_limits<double>::quiet_NaN();

  probs.resize(mod.fpmm_sub_model.loci());

  for(size_t l = 0; l < mod.fpmm_sub_model.loci(); ++l) {
    probs[l].resize((l+1)*2 + 1);
    for(size_t x = 0; x < ((l+1)*2)+1; ++x) { 
      probs[l][x].resize((l+1)*2 + 1);
      for(size_t y = 0; y < ((l+1)*2)+1; ++y) {
        probs[l][x][y].resize((l+1)*2 + 1);
        for(size_t z = 0; z < ((l+1)*2)+1; ++z) {
          probs[l][x][y][z] = QNAN;
        }
      }
    }
  }

  // Calculate probs

  for(size_t l = 0; l < mod.fpmm_sub_model.loci(); ++l) {
    for(size_t x = 0; x < (l+1)*2+1; ++x) { 
      for(size_t y = 0; y < (l+1)*2+1; ++y) {
        for(size_t z = 0; z < (l+1)*2+1; ++z) {
          probs[l][x][y][z] = int_prob(x,y,z,l+1);
/*
          cout << "[l=" << l
               << "][Vi=" << x
               << "][Vm=" << y
               << "][Vf=" << z
               << "]= " << probs[l][x][y][z] << endl;
*/
        }
      }
    }
  }

  // Release un-needed memory:
/*
  for(size_t l = 0; l < mod.fpmm_sub_model.loci() - 1; ++l)
    probs[l].clear();
*/
}

//==========================================================================
// prob(...)
//==========================================================================
double
polygenic_transition_calculator::prob(size_t polygenotype_indiv,
                                      size_t polygenotype_mother, 
                                      size_t polygenotype_father,
    				      size_t num_of_locii)
{
  return probs[num_of_locii-1     ][polygenotype_indiv ]
              [polygenotype_mother][polygenotype_father];
}

//==========================================================================
// int_prob(...)
//==========================================================================
double
polygenic_transition_calculator::int_prob(size_t polygenotype_indiv,
                                          size_t polygenotype_mother, 
                                          size_t polygenotype_father,
    				          size_t num_of_locii)
{
  size_t Vi = polygenotype_indiv;
  size_t Vm = polygenotype_mother;
  size_t Vf = polygenotype_father;

  if(polygenotype_indiv > num_of_locii * 2) return 0;
   else Vi = polygenotype_indiv;

  if(polygenotype_mother > num_of_locii * 2) return 0;
   else Vm = polygenotype_mother;

  if(polygenotype_father > num_of_locii * 2) return 0;
   else Vf = polygenotype_father;

  if(SAGE::isnan(probs[num_of_locii-1][Vi][Vm][Vf]))
    probs[num_of_locii-1][Vi][Vm][Vf] = calculate_prob(Vi,Vm,Vf,num_of_locii);

  return probs[num_of_locii-1][Vi][Vm][Vf];
}

//==========================================================================
// calculate_prob(...)
//==========================================================================
double
polygenic_transition_calculator::calculate_prob(size_t polygenotype_indiv,
                                                size_t polygenotype_mother, 
                                                size_t polygenotype_father,
                                                size_t num_of_locii)
{
  static const double taus[] = { 0.0, 0.5, 1.0 };

  double polygenic_prob = 0;

  if(num_of_locii == 1)
  {
    genotype_index genotype_indiv  = (genotype_index)polygenotype_indiv;
    genotype_index genotype_mother = (genotype_index)polygenotype_mother;
    genotype_index genotype_father = (genotype_index)polygenotype_father;

    double tm = taus[genotype_mother];
    double tf = taus[genotype_father];

    switch(genotype_indiv)
    {
      case index_AA : polygenic_prob = tm * tf;                           break;
      case index_AB : polygenic_prob = tm * (1.0 - tf) + tf * (1.0 - tm); break;
      case index_BB : polygenic_prob = (1.0 - tm) * (1.0 - tf);           break;
    }
  }
  else
  {
    double trans1                     = 0.0;
    double trans2                     = 0.0;
    double trans3                     = 0.0;
    double trans4                     = 0.0;
    double genotype_indiv_loop_total  = 0.0;
    double genotype_mother_loop_total = 0.0;
    double genotype_father_loop_total = 0.0;
    size_t new_polygenotype_indiv     = 0;
    size_t new_polygenotype_mother    = 0;
    size_t new_polygenotype_father    = 0;

    genotype_indiv_loop_total = 0;
    for(size_t genotype_indiv_loop = 0; 
              (genotype_indiv_loop < 3) && (polygenotype_indiv >= genotype_indiv_loop); 
             ++genotype_indiv_loop)
    {
      genotype_mother_loop_total = 0;
      for(size_t genotype_mother_loop = 0; 
                (genotype_mother_loop < 3) && (polygenotype_mother >= genotype_mother_loop); 
               ++genotype_mother_loop)
      {
        genotype_father_loop_total = 0;
        for(size_t genotype_father_loop = 0; 
                  (genotype_father_loop < 3) && (polygenotype_father >= genotype_father_loop); 
                 ++genotype_father_loop)
        { 
          new_polygenotype_indiv  = polygenotype_indiv  - genotype_indiv_loop;
          new_polygenotype_mother = polygenotype_mother - genotype_mother_loop;
          new_polygenotype_father = polygenotype_father - genotype_father_loop;  

          // These four elements are representative of the four elements in
          // the equation near equation 80 in the Formula document.  trans1
          // is the first, trans2 is the second, etc.

          trans1 = int_prob(new_polygenotype_indiv,new_polygenotype_mother,
                   new_polygenotype_father,num_of_locii-1);
          trans2 = int_prob(genotype_indiv_loop, genotype_mother_loop, 
                   genotype_father_loop, 1);
          trans3 = bin_trans(genotype_mother_loop,new_polygenotype_mother,polygenotype_mother,num_of_locii);
          trans4 = bin_trans(genotype_father_loop,new_polygenotype_father,polygenotype_father,num_of_locii);

          genotype_father_loop_total += trans1 * trans2 * trans3 * trans4;        

        } // End of genotype father loop

        genotype_mother_loop_total += genotype_father_loop_total;

      } // End of genotype mother loop

      genotype_indiv_loop_total += genotype_mother_loop_total;    

    } // End of genotype indiv loop

    polygenic_prob = genotype_indiv_loop_total;
  }

  if(SAGE::isnan(polygenic_prob)) exit(0);

  return polygenic_prob;
}
   
//========================================================================== 
// bin_trans(...)
//==========================================================================
double
polygenic_transition_calculator::bin_trans(size_t genotype,  
	size_t new_polygenotype, size_t old_polygenotype, size_t num_of_locii)
{
  double bin_trans;

  double p = mod.fpmm_sub_model.frequency(); // See eq 80 for explanation of this default value

  // Calculate numerator:

  double bin1 = NUMERICS::bin_prob(2*(num_of_locii-1),p,new_polygenotype);
  double bin2 = NUMERICS::bin_prob(2,p,genotype);                         

  // Calculate denominator:

  size_t i;
  double denominator = 0;
  double bin3;
  double bin4;

  for(size_t j = 0; j < 3; ++j)
  {
    i            = old_polygenotype - j;
    bin3         = NUMERICS::bin_prob(2 * (num_of_locii-1),p,i);  
    bin4         = NUMERICS::bin_prob(2,p,j);                     
    denominator += bin3 * bin4;
  }

  bin_trans = (bin1 * bin2) / denominator;
  return bin_trans;
}

//========================================================================== 
// get_num_of_locii() 
//==========================================================================
size_t
polygenic_transition_calculator::get_num_of_locii()
{ return 5; }

}
}
