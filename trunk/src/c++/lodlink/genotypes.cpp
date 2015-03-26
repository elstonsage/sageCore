//============================================================================
// File:      genotypes.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   5/21/3 - created.                         djb
//                                                                          
// Notes:     Implementation genotype probability classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/genotypes.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  genotype_probs
//============================================================================
//
void
genotype_probs::calculate_genotypes(peeler& plr, const member_type& ind, genotype_result& result)
{
  MLOCUS::penetrance_model  trait_pm(ge_models::get_model(plr.subpedigree(), plr.trait()));
  MLOCUS::penetrance_model  marker_pm(ge_models::get_model(plr.subpedigree(), plr.marker()));
  size_t  trait_phenotype(ind.subindex() + 1);
  size_t  marker_phenotype(ind.subindex() + 1);
  
  log_double  likes[] = { log_double(0), log_double(0), log_double(0) };
  log_double  total_like(0);
  
  MLOCUS::penetrance_model::phased_penetrance_iterator  t_iter = trait_pm.phased_penetrance_begin(trait_phenotype);
  for(; t_iter != trait_pm.phased_penetrance_end(trait_phenotype); ++t_iter)
  {
    log_double  like(0);
  
    MLOCUS::penetrance_model::phased_penetrance_iterator  m_iter = marker_pm.phased_penetrance_begin(marker_phenotype);
    for(; m_iter != marker_pm.phased_penetrance_end(marker_phenotype); ++m_iter)
    {
      log_double like_term(1);
      like_term *= plr.anterior(ind, joint_pen_iter(t_iter, m_iter));
      like_term *= joint_pen_iter(t_iter, m_iter).penetrance();
      like_term *= plr.posterior(ind, joint_pen_iter(t_iter, m_iter));  
    
      like += like_term;
    }
    
    total_like += like;
    
    string  genotype = t_iter.phased_geno().name();
    if(genotype == "A<A")
    {
      likes[aa] += like;
    }
    else if(genotype == "A<B")
    {
      likes[ab] += like;
    }
    else if(genotype == "B<A")
    {
      likes[ab] += like;
    } 
    else if(genotype == "B<B")
    {
      likes[bb] += like;
    }
    else
    { 
      assert(false);
    }       
  }
  
  geno_probs  temp_probs;
  
  temp_probs.aa = (likes[aa] / total_like).get_double();
  temp_probs.ab = (likes[ab] / total_like).get_double();
  temp_probs.bb = (likes[bb] / total_like).get_double();
  
  assert(abs(temp_probs.aa + temp_probs.ab + temp_probs.bb - 1) < .0001);
  
  assert(result.probs.insert(genotype_result::geno_map::value_type(&ind, temp_probs)).second);
}

void
genotype_probs::write_summary(ostream& out) const
{
  // INTENTIONALLY LEFT BLANK
}


//============================================================================
// IMPLEMENTATION:  non_ss_genotype_probs
//============================================================================
//
void
non_ss_genotype_probs::calculate()
{
  do_task_calculations<non_ss_genotype_probs, non_ss_genotype_result>(*this);
}

// - Here the name, calculate_alt(), is a misnomer, but it is expedient  
//   to use do_task_calculations(), which calls a function by this name.
//
// - Estimate the recombination fraction between the trait and the current marker 
//   using all data in the pedigree data file.
//   For every extended family (subpedigree).
//      For every individual in the extended family.
//         For each of 3 type genotypes for the current individual
//            Calculate the likelihood for the extended family
//         Store the likelihood of each genotype divided by the sum of all three likelihoods.
//   
void
non_ss_genotype_probs::calculate_alt(size_t trait_index, size_t marker_index, non_ss_genotype_result& result)
{
  result.theta = estimate_theta(trait_index, marker_index);
  
  if(SAGE::isnan(result.theta))
  {
    return;
  }

  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      mle_sub_model  mle(false, false);
      mle.set_strict_limits();
      mle.set_average_theta(result.theta);
      
      peeler  plr(*subped_iter, mle, trait_index, marker_index);
      member_const_iterator  member_iter = subped_iter->member_begin();
      for(; member_iter != subped_iter->member_end(); ++member_iter)
      {
        calculate_genotypes(plr, *member_iter, result);
      }      
    }
  }
}

double  
non_ss_genotype_probs::estimate_theta(size_t trait_index, size_t marker_index)
{
  mle_sub_model  mle(false, false);
  mle.set_strict_limits();
  
  if(likelihood_finite(mle, trait_index, marker_index, "genotype probability calculations"))
  {
    MAXFUN::Results  data = maximize(mle, trait_index, marker_index);
    if(! SAGE::isnan(max_checked_value(data)))
    {
      return  mle.average_theta();
    }
  }
  
  return  QNAN;
}

void
non_ss_genotype_probs::calculate_null(size_t trait_index, size_t marker_index, non_ss_genotype_result& result)
{
  // INTENTIONALLY LEFT BLANK
}

void
non_ss_genotype_probs::write_detail(ostream& out) const
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
      member_const_iterator  member_iter;
      for(member_iter = subped_iter->member_begin();
          member_iter != subped_iter->member_end();
          ++member_iter)
      {
        out << setw(non_ss_genotype_result::member_header_offset) << ""
            << "Pedigree " << ped_iter->name() 
            << " Member " << member_iter->name() << endl;
            
        out << non_ss_genotype_result::detail_columns;
        for(size_t i = 0; i < my_results.size(); ++i)
        {
          non_ss_genotype_result*  ptr = static_cast<non_ss_genotype_result*>(my_results[i].get());
          ptr->write_genotypes(out, *member_iter);
        }
        
        out << "\n" << endl;
      }
    }
  }  
}

void
non_ss_genotype_probs::write_detail_header(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << "\n\n" << endl;
  out << setfill(' ');
  out << left;
    
  out << setw(_LOCUS.tw()) << ""
      << GENOTYPES << " " << "\n\n" << endl;
  
  out.flags(old_flags); 
}


//============================================================================
// IMPLEMENTATION:  ss_genotype_probs
//============================================================================
//
void
ss_genotype_probs::calculate()
{
  do_task_calculations<ss_genotype_probs, ss_genotype_result>(*this);
}

// - Here the name, calculate_alt(), is a misnomer, but it is expedient  
//   to use do_task_calculations(), which calls a function by this name.
//
// - Estimate the recombination fractions between the trait and the current marker 
//   using all data in the pedigree data file.
//   For every extended family (subpedigree).
//      For every individual in the extended family.
//         For each of 3 type genotypes for the current individual
//            Calculate the likelihood for the extended family
//         Store the likelihood of each genotype divided by the sum of all three likelihoods.
//   
void
ss_genotype_probs::calculate_alt(size_t trait_index, size_t marker_index, ss_genotype_result& result)
{
  result.thetas = estimate_thetas(trait_index, marker_index);
  
  if(SAGE::isnan(result.thetas.male_theta) || SAGE::isnan(result.thetas.female_theta))
  {
    return;
  }

  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
  
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      mle_sub_model  mle(true, false);
      mle.set_strict_limits();
      mle.set_male_theta(result.thetas.male_theta);
      mle.set_female_theta(result.thetas.female_theta);      
      
      peeler  plr(*subped_iter, mle, trait_index, marker_index);

      member_const_iterator  member_iter = subped_iter->member_begin();
      for(; member_iter != subped_iter->member_end(); ++member_iter)
      {
        calculate_genotypes(plr, *member_iter, result);
      }      
    }
  }
}

theta_pair  
ss_genotype_probs::estimate_thetas(size_t trait_index, size_t marker_index)
{
  mle_sub_model  mle(true, false);
  mle.set_strict_limits();
  
  if(likelihood_finite(mle, trait_index, marker_index, "genotype probability calculations"))
  {
    MAXFUN::Results  data = maximize(mle, trait_index, marker_index);
    if(! SAGE::isnan(max_checked_value(data)))
    {
      return  theta_pair(mle.male_theta(), mle.female_theta());
    }
  }
  
  return  theta_pair(QNAN, QNAN);
}

void
ss_genotype_probs::calculate_null(size_t trait_index, size_t marker_index, ss_genotype_result& result)
{
  // INTENTIONALLY LEFT BLANK
}

void
ss_genotype_probs::write_detail(ostream& out) const
{
  write_detail_header(out);

  for(pedigree_const_iterator ped_iter = my_mped.pedigree_begin(); ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    for(subpedigree_const_iterator subped_iter = ped_iter->subpedigree_begin(); subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      for(member_const_iterator member_iter = subped_iter->member_begin(); member_iter != subped_iter->member_end(); ++member_iter)
      {
        out << setw(ss_genotype_result::member_header_offset) << ""
            << "Pedigree " << ped_iter->name() 
            << " Member " << member_iter->name() << endl;
            
        out << ss_genotype_result::detail_columns;
        for(size_t i = 0; i < my_results.size(); ++i)
        {
          ss_genotype_result*  ptr = static_cast<ss_genotype_result*>(my_results[i].get());
          ptr->write_genotypes(out, *member_iter);
        }
        
        out << "\n" << endl;
      }
    }
  }  
}

void
ss_genotype_probs::write_detail_header(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << "\n\n" << endl;
  out << setfill(' ');
  out << left;
    
  out << setw(_LOCUS.tw()) << ""
      << GENOTYPES << " " << "\n\n" << endl;
  
  out.flags(old_flags); 
}

}
}
