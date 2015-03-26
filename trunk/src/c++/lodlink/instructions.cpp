//============================================================================
// File:      instructions.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   9/12/2 - created.                         djb
//                                                                          
// Notes:     Default lods.  Implementaion of instructions class.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/instructions.h"

namespace SAGE
{

namespace LODLINK
{

const double  THETAS[] = { 0, .01, .05, .1, .2, .3, .4 };
const theta_pair  THETA_PAIRS[] = { theta_pair(0, 0),     theta_pair(.01, .01),
                                    theta_pair(.05, .05), theta_pair(.1, .1),
                                    theta_pair(.2, .2),   theta_pair(.3, .3),
                                    theta_pair(.4, .4)                         };

const double*  FIRST_THETA = THETAS;
const double*  LAST_THETA  = THETAS + sizeof(THETAS) / sizeof(THETAS[0]);
const theta_pair*  FIRST_THETA_PAIR = THETA_PAIRS;
const theta_pair*  LAST_THETA_PAIR  = THETA_PAIRS + sizeof(THETA_PAIRS) / sizeof(THETA_PAIRS[0]);

//============================================================================
// IMPLEMENTATION:  instructions
//============================================================================
//
void
instructions::write(ostream& out) const
{
  bool  lod_scores = ! (average_thetas.empty() && male_female_thetas.empty());

  ios::fmtflags old_flags = out.flags();
  out << left;
  out << setfill(' ');
  
  out << "Options Selected\n"
      << "================" << endl;
      
  out << setw(OPTION_SZ) << "Main locus type"
      << link_option() << "\n"
      << setw(OPTION_SZ) << "Main locus name"
      << trait << "\n\n"
      
      << setw(OPTION_SZ) << "Lod scores"
      << (lod_scores ? "yes" : "no") << "\n\n";
      
 if(lod_scores)
 {
   write_thetas(out);
 }
      
  out << setw(OPTION_SZ) << "Linkage tests"
      << (linkage_test ? "yes" : "no") << "\n";
      
  if(linkage_test)
  {
    out << setw(OPTION_SZ) << "  Sex-specific recombination fractions"
        << (linkage_sex_specific ? "yes" : "no") << "\n"
        << setw(OPTION_SZ) << "  Assume homogeneity"
        << (linkage_homog ? "yes" : "no") << "\n\n";
  }
  else
  {
    out << "\n";
  }
     
  out << setw(OPTION_SZ) << "Smith's test for homogeneity"
      << (smiths_test ? "yes" : "no") << "\n";
      
  if(smiths_test)
  {
    out << setw(OPTION_SZ) << "  Sex-specific recombination fractions"
        << (smiths_sex_specific ? "yes" : "no") << "\n\n";
  }
  else
  {
    out << "\n";
  }  
      
  out << setw(OPTION_SZ) << "Morton's test for homogeneity"
      << (mortons_test ? "yes" : "no") << "\n\n";
      
  if(mortons_test)
  {
    out << setw(OPTION_SZ) << "  Sex-specific recombination fractions"
        << (mortons_sex_specific ? "yes" : "no") << "\n\n";
    write_groups(out);
  }
       
  out << setw(OPTION_SZ) << "Genotype probabilities"
      << (genotypes ? "yes" : "no") << "\n";
      
  if(genotypes)
  {
    out << setw(OPTION_SZ) << "  Sex-specific recombination fractions"
        << (genotypes_sex_specific ? "yes" : "no") << "\n\n" << endl;      
  }
  else
  {
    out << "\n" << endl;
  }
      
  out.flags(old_flags);
}

void
instructions::write_groups(ostream& out) const
{
  const size_t  LABEL_SZ = 35;

  ios::fmtflags old_flags = out.flags();
  out << left << setfill(' '); 
  out << std::fixed << setprecision(PRC1);

  out << "  Groups Selected\n"
      << "  ===============" << endl;
    
  map<string, group>::const_iterator  g_iter;
  for(g_iter = groups.begin(); g_iter != groups.end(); ++g_iter)
  {
    string  label = "  group " + g_iter->first + " pedigrees";
    out << setw(LABEL_SZ) << label;
    group::const_iterator  p_iter;
    for(p_iter = g_iter->second.begin(); p_iter != g_iter->second.end(); ++p_iter)
    {
      out << *p_iter << setw(SPACE_SZ) << ""; 
    }
    
    out << endl;
  }
    
  out << endl;
  
  out.flags(old_flags);
}

void
instructions::write_thetas(ostream& out) const
{
  const size_t  LABEL_SZ = 15;

  ios::fmtflags old_flags = out.flags();
  out << left << setfill(' '); 
  out << std::fixed << setprecision(PRC1);

  out << "  Recombination Fractions Selected\n"
      << "  ================================" << endl;

  if(! average_thetas.empty())
  {
    assert(male_female_thetas.empty());
    
    out << setw(LABEL_SZ) << "  Sex averaged";
    
    out << right << "  ";
    for(size_t i = 0; i < average_thetas.size(); ++i)
    {
      out << setw(RECOM_SZ) << average_thetas[i]
          << setw(SPACE_SZ) << "";    
    }
    
    out << "\n" << endl;
  }
  
  if(! male_female_thetas.empty())
  {
    assert(average_thetas.empty());
    
    out << setw(LABEL_SZ) << "  Male";
    
    out << right << "  ";
    for(size_t i = 0; i < male_female_thetas.size(); ++i)
    {
      out << setw(RECOM_SZ) << male_female_thetas[i].male_theta
          << setw(SPACE_SZ) << "";    
    }
    
    out << "\n"
        << left << setw(LABEL_SZ) << "  Female";
    
    out << right << "  ";
    for(size_t i = 0; i < male_female_thetas.size(); ++i)
    {
      out << setw(RECOM_SZ) << male_female_thetas[i].female_theta
          << setw(SPACE_SZ) << "";    
    }
    
    out << "\n" << endl;
  }
  
  out.flags(old_flags);
}

}
}

