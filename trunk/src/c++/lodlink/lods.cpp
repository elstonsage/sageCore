//============================================================================
// File:      lods.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   2/13/3 - created.                         djb
//                                                                          
// Notes:     Implementation of lod score calculation classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/lods.h"

namespace SAGE
{

namespace LODLINK
{


//============================================================================
// IMPLEMENTATION:  non_ss_lods
//============================================================================
//
void
non_ss_lods::calculate()
{
  do_task_calculations<non_ss_lods, non_ss_lods_result>(*this);
}

// - Here the name, calculate_alt(), is a misnomer, but it is expedient  
//   to use do_task_calculations(), which calls a function by this name.
//
void
non_ss_lods::calculate_alt(size_t trait_index, size_t marker_index, non_ss_lods_result& result)
{
  mle_sub_model  mle(false, false);
  for(size_t i = 0; i < my_instructions.average_thetas.size(); ++i)
  {
    double  theta = my_instructions.average_thetas[i];
    mle.set_average_theta(theta);

    multipedigree_lod_score  mped_struct;
    double  mped_ls = 0;

    pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
    for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
    {
      string  ped_name = ped_iter->name();
      pedigree_lod_score  ped_ls_struct(ped_name);
      double  ped_ls = 0;
    
      subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
      for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
      {
        peeler  p(*subped_iter, mle, trait_index, marker_index);
        subped_calculator  sp_calc(p);
        double  alt_ln_like  = sp_calc.likelihood().get_log();  
        double  null_ln_like = sp_calc.unlinked_likelihood().get_log(); 
        double  subped_ls(lod_score(alt_ln_like, null_ln_like)); 
        
        string  member_name = subped_iter->member_begin()->name();
        
        ped_ls += subped_ls;
        ped_ls_struct.sub_lod_scores.push_back(subpedigree_lod_score(member_name, subped_ls)); 
      }
      
      ped_ls_struct.lod_score = ped_ls;
      mped_struct.ped_lod_scores.push_back(ped_ls_struct);
      mped_ls += ped_ls;
    }
    
    mped_struct.lod_score = mped_ls;
    result.lod_scores.push_back(make_pair(theta, mped_struct));
  }
}

// - Doesn't need to do anything.
//
void
non_ss_lods::calculate_null(size_t trait_index, size_t marker_index, non_ss_lods_result& result)
{}

void
non_ss_lods::write_summary(ostream& out) const
{
  write_summary_headers(out);
  
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
non_ss_lods::write_summary_headers(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << left;
  for(int i = 0; i < 3; ++i)
  {
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << setw(_LODS_AVE_.lw()) << _LODS_AVE_[i] << "\n";
  }    
  
  out << endl;
  write_columns_header(out);
  
  out.flags(old_flags); 
}

void
non_ss_lods::write_columns_header(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();

  out << setfill(' ');
  out << left;

  out << setw(_LOCUS.tw()) << _LOCUS[0];
  out << std::fixed << setprecision(PRC1);
  for(size_t i = 0; i < my_instructions.average_thetas.size(); ++i)
  {
    out << setw(RECOM_SZ) << my_instructions.average_thetas[i]
        << setw(SPACE_SZ) << "";
  }
  
  out << endl;
  out << setfill(U_CHR) << setw(_LOCUS.lw()) << "";
  out << setfill(' ') << setw(_LOCUS.sw()) << "";
  for(size_t i = 0; i < my_instructions.average_thetas.size(); ++i)
  {
    out << setfill(U_CHR) << setw(RECOM_SZ) << ""
        << setfill(' ') << setw(SPACE_SZ) << "";
  }
  
  out << endl;
  
  out.flags(old_flags);
}

void
non_ss_lods::write_detail(ostream& out) const
{
  write_detail_header(out);

  pedigree_const_iterator  ped_iter;
  for(ped_iter = my_mped.pedigree_begin(); ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    subpedigree_const_iterator  subped_iter;
    for(subped_iter = ped_iter->subpedigree_begin(); 
        subped_iter != ped_iter->subpedigree_end();
        ++subped_iter)
    {
      out << setw(_LOCUS.tw()) << ""
          << "Constituent Pedigree in Pedigree " << ped_iter->name() 
          << " Containing Member " << subped_iter->member_begin()->name()
          << "\n" << endl;
          
      write_columns_header(out);
      for(size_t i = 0; i < my_results.size(); ++i)
      {
        non_ss_lods_result*  ptr = static_cast<non_ss_lods_result*>(my_results[i].get());
        ptr->write_family_detail(out, ped_iter->name(), subped_iter->member_begin()->name());
      }
      
      out << "\n" << endl;
    }
  }
}

void
non_ss_lods::write_detail_header(ostream& out)
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << left;
  for(int i = 0; i < 3; ++i)
  {
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << setw(_LODS_AVE_FAM_.lw()) << _LODS_AVE_FAM_[i] << "\n";
  }    
  
  out << endl;
  
  out.flags(old_flags); 
}

//============================================================================
// IMPLEMENTATION:  ss_lods
//============================================================================
//
void
ss_lods::calculate()
{
  do_task_calculations<ss_lods, ss_lods_result>(*this);
}

// - Here the name, calculate_alt(), is a misnomer, but it is expedient  
//   to use do_task_calculations(), which uses this name.
//
void
ss_lods::calculate_alt(size_t trait_index, size_t marker_index, ss_lods_result& result)
{
  mle_sub_model  mle(true, false);
  for(size_t i = 0; i < my_instructions.male_female_thetas.size(); ++i)
  {
    double  male_theta = my_instructions.male_female_thetas[i].male_theta;
    double  female_theta = my_instructions.male_female_thetas[i].female_theta;
    mle.set_male_theta(male_theta);
    mle.set_female_theta(female_theta);
    
    multipedigree_lod_score  mped_struct;

    double  mped_ls = 0;

    pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
    for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
    {
      string  ped_name = ped_iter->name();
      pedigree_lod_score  ped_ls_struct(ped_name);
      double  ped_ls = 0;
    
      subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
      for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
      {
        peeler  p(*subped_iter, mle, trait_index, marker_index);
        subped_calculator  sp_calc(p);
        double  alt_ln_like  = sp_calc.likelihood().get_log();
        double  null_ln_like = sp_calc.unlinked_likelihood().get_log();
        double  subped_ls(lod_score(alt_ln_like, null_ln_like));
        
        string  member_name = subped_iter->member_begin()->name();
        
        ped_ls += subped_ls;
        ped_ls_struct.sub_lod_scores.push_back(subpedigree_lod_score(member_name, subped_ls)); 
      }
      
      ped_ls_struct.lod_score = ped_ls;
      mped_struct.ped_lod_scores.push_back(ped_ls_struct);
      mped_ls += ped_ls;
    }
    
    mped_struct.lod_score = mped_ls;
    result.lod_scores.push_back(make_pair(theta_pair(male_theta, female_theta), mped_struct));
  }
}

// - Doesn't need to do anything.
//
void
ss_lods::calculate_null(size_t trait_index, size_t marker_index, ss_lods_result& result)
{}

void
ss_lods::write_summary(ostream& out) const
{
  write_summary_headers(out);
  
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
ss_lods::write_summary_headers(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << left;
  for(int i = 0; i < 3; ++i)
  {
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << setw(_LODS_SEX_SPEC_.lw()) << _LODS_SEX_SPEC_[i] << "\n";
  }    
  
  out << endl;
  write_columns_header(out);
  
  out.flags(old_flags);
}

void
ss_lods::write_columns_header(ostream& out) const
{
  ios::fmtflags  old_flags = out.flags();
  
  out << setfill(' ');
  out << left;
  
  out << std::fixed << setprecision(PRC1);
  
  out << setw(_LOCUS.tw()) << "";
  for(size_t i = 0; i < my_instructions.male_female_thetas.size(); ++i)
  {
    out << setw(RECOM_SZ) << my_instructions.male_female_thetas[i].male_theta
        << setw(SPACE_SZ) << "";
  }
  
  out << endl;
  
  out << setw(_LOCUS.tw()) << _LOCUS[0];
  for(size_t i = 0; i < my_instructions.male_female_thetas.size(); ++i)
  {
    out << setw(RECOM_SZ) << my_instructions.male_female_thetas[i].female_theta
        << setw(SPACE_SZ) << "";
  }
  
  out << endl;
  
  out << setfill(U_CHR) << setw(_LOCUS.lw()) << "";
  out << setfill(' ') << setw(_LOCUS.sw()) << "";
  for(size_t i = 0; i < my_instructions.male_female_thetas.size(); ++i)
  {
    out << setfill(U_CHR) << setw(RECOM_SZ) << ""
        << setfill(' ') << setw(SPACE_SZ) << "";
  }
  
  out << endl;
  
  out.flags(old_flags);
}

void
ss_lods::write_detail(ostream& out) const
{
  write_detail_header(out);

  pedigree_const_iterator  ped_iter;
  for(ped_iter = my_mped.pedigree_begin(); ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    subpedigree_const_iterator  subped_iter;
    for(subped_iter = ped_iter->subpedigree_begin(); 
        subped_iter != ped_iter->subpedigree_end();
        ++subped_iter)
    {
      out << setw(_LOCUS.tw()) << ""
          << "Constituent Pedigree in Pedigree " << ped_iter->name() 
          << " Containing Member " << subped_iter->member_begin()->name()
          << "\n" << endl;
          
      write_columns_header(out);
      for(size_t i = 0; i < my_results.size(); ++i)
      {
        ss_lods_result*  ptr = static_cast<ss_lods_result*>(my_results[i].get());
        ptr->write_family_detail(out, ped_iter->name(), subped_iter->member_begin()->name());
      }
      
      out << "\n" << endl;
    }
  }
}

void
ss_lods::write_detail_header(ostream& out)
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ');
  out << left;
  for(int i = 0; i < 3; ++i)
  {
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << setw(_LODS_SEX_SPEC_FAM_.lw()) << _LODS_SEX_SPEC_FAM_[i] << "\n";
  }    
  
  out << endl;
  
  out.flags(old_flags); 
}

}
}
