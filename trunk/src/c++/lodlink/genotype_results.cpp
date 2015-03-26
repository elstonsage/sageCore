//============================================================================
// File:      genotype_results.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   5/21/3 - created.                         djb
//                                                                          
// Notes:     Implementation of genotype probabilities results classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/genotype_results.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  non_ss_genotype_result
//============================================================================
//
size_t  non_ss_genotype_result::member_header_offset;
header  non_ss_genotype_result::detail_columns;

void
non_ss_genotype_result::build_headers()
{
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_THREE_LINE_LOCUS);
  detail_columns.add_col(_RECOM_SMALL_);  
  detail_columns.add_col(_TYPE_AA);  
  detail_columns.add_col(_TYPE_AB);  
  detail_columns.add_col(_TYPE_BB);
  
  member_header_offset = detail_columns[0].tw() +
                         detail_columns[1].tw()  ;
}

void  
non_ss_genotype_result::write_genotypes(ostream& out, const member& ind) const
{
  map_iter  iter = probs.find(&ind);
  assert(iter != probs.end());
  
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(non_ss_genotype_result::detail_columns.col_w(1)); 
      write_double(out, theta); 
  out << setw(non_ss_genotype_result::detail_columns.spc_w(1)) << "";      
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);

  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.aa);
  out << setw(non_ss_genotype_result::detail_columns.spc_w(2)) << ""; 
  
  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.ab);
  out << setw(non_ss_genotype_result::detail_columns.spc_w(3)) << "";
  
  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.bb);
  out << setw(non_ss_genotype_result::detail_columns.spc_w(4)) << "";      
  
  out << endl;
  
  out.flags(old_flags);  
}


//============================================================================
// IMPLEMENTATION:  ss_genotype_result
//============================================================================
//
size_t  ss_genotype_result::member_header_offset;
header  ss_genotype_result::detail_columns;

void
ss_genotype_result::build_headers()
{
  detail_columns.set_underline(U_CHR);
  detail_columns.add_col(_THREE_LINE_LOCUS);  
  detail_columns.add_col(_M_RECOM_SMALL_);
  detail_columns.add_col(_F_RECOM_SMALL_);  
  detail_columns.add_col(_TYPE_AA);  
  detail_columns.add_col(_TYPE_AB);  
  detail_columns.add_col(_TYPE_BB);  
  
  member_header_offset = detail_columns[0].tw() +
                         detail_columns[1].tw() +
                         detail_columns[2].tw()  ;
}

void  
ss_genotype_result::write_genotypes(ostream& out, const member& ind) const
{
  map_iter  iter = probs.find(&ind);
  assert(iter != probs.end());
  
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << ""
      
      << right << std::fixed << setprecision(PRC1)
      << setw(ss_genotype_result::detail_columns.col_w(1)); 
      write_double(out, thetas.male_theta); 
  out << setw(ss_genotype_result::detail_columns.spc_w(1)) << ""      
  
      << setw(ss_genotype_result::detail_columns.col_w(2)); 
      write_double(out, thetas.female_theta); 
  out << setw(ss_genotype_result::detail_columns.spc_w(2)) << "";
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);

  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.aa);
  out << setw(ss_genotype_result::detail_columns.spc_w(3)) << ""; 
  
  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.ab);
  out << setw(ss_genotype_result::detail_columns.spc_w(4)) << "";
  
  out << setw(GENOTYPE_SZ);
  write_double(out, iter->second.bb);
  out << setw(ss_genotype_result::detail_columns.spc_w(5)) << "";      
  
  out << endl;
  
  out.flags(old_flags);  
}

}
}
