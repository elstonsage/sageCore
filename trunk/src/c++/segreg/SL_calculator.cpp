//=============================================================
//
//  File:	SL_calculator.cpp
//
//  Authors:	Stephen Gross
//		Kai He
//
//  History:	sag/kh	Initial implementation	Jul 23 01
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//
//=============================================================

#include "segreg/SL_calculator.h"

using namespace std;
namespace SAGE {
namespace SEGREG {

//============================================================
//
// FPMM_SL::FPMM_SL(...)
//
//============================================================
FPMM_SL::FPMM_SL(
    const FPED::Multipedigree & ped_data,
    const model               & modp,
    bool                        use_ascertainmentp) :
    trans_poly_calc(modp),
    pen(ped_data,modp,use_ascertainmentp),
    mod(modp)
{ }

//=============================================================
//
//  FPMM_SL::int_founder_SL2_4(...)
//
//=============================================================
log_double
FPMM_SL::int_founder_SL2_4(const penetrance_info&	    indiv) const
{
  // SL2(i, P) = psi(Ui,Vi)P(Ti|Ui,Vi)

  log_double SL2(1.0);

  double psi        = mod.freq_sub_model.prob    (indiv.genotype) * 
                      mod.fpmm_sub_model.pop_freq(indiv.polygenotype);
  double penetrance = pen.get_polygenic_penetrance(indiv);

  SL2 = psi * penetrance;

  return SL2;
}

//============================================================
//
//  FPMM_SL::int_founder_SL6(...)
//
//============================================================
log_double
FPMM_SL::int_founder_SL6(const penetrance_info&	   indiv,
		         const penetrance_info&	   spouse) const
{
  // SL6(i, s, P) is to be calculated the same way as for a non-founder.
  // See equation #26

  log_double SL6 = nonfounder_SL6(indiv,spouse);
  return SL6;
}

/**@file
  * \f[ SL2(i,P)=
  *   P(u_i|u_m,u_f) P(v_i|v_m,v_f) P(t_i|u_i,v_i) \prod_{b \in B_i}
  *     \left\{
  *       \sum_{u_b} \sum_{v_b} 
  *       P(u_b|u_m,u_f) P(v_b|v_m,v_f) P(t_b|u_b,v_b)
  *       \prod_{h \in S_b} pos_bh(u_b,v_b)
  *     \right\}
  * \f]
  *
  *
  *
  * \f[ SL4(i,s,P)=
  *   P(u_i|u_m,u_f) P(v_i|v_m,v_f) P(t_i|u_i,v_i,u_s,v_s,t_s) 
  *   \prod_{b \in B_i}
  *     \left\{
  *       \sum_{u_b} \sum_{v_b} 
  *       P(u_b|u_m,u_f) P(v_b|v_m,v_f) P(t_b|u_b,v_b)
  *       \prod_{h \in S_b} pos_bh(u_b,v_b)
  *     \right\}
  * \f]
  *
  *
  *
  * \f[ SL6(i,s,P)=
  *   \prod_{k \in C_is}
  *     \left\{
  *       \sum_{u_k} \sum_{v_k} 
  *       P(u_k|u_i,u_s) P(v_k|v_i,v_s) P(t_k|u_k,v_k)
  *       \prod_{h \in S_k} pos_kh(u_k,v_k)
  *     \right\}
  * \f]
  */

//===================================================================
//
//  FPMM_SL::int_nonfounder_SL2_4(...)
//
//===================================================================

log_double
FPMM_SL::int_nonfounder_SL2_4(const penetrance_info&          indiv,
		              const penetrance_info&          mother,
		              const penetrance_info&          father) const
{
  // For definition of SL2 and SL4, see equations #38, 39.

  log_double SL2_4(1.0);

  penetrance_info        sibling;

  // Transmissions (P(u_i | u_m, u_f) and P(v_i | v_m, v_f))

  double trans_geno_indiv     = mod.transm_sub_model.prob (indiv .genotype,
                                                           mother.genotype,
                                                           father.genotype);
  double trans_polygeno_indiv = trans_poly_calc     .prob
           (indiv .polygenotype, mother.polygenotype, father.polygenotype,
            mod.fpmm_sub_model.loci());

  if(!trans_geno_indiv || !trans_polygeno_indiv)
    return log_double(0.0);

  SL2_4 *= trans_geno_indiv * trans_polygeno_indiv;

  // Penetrance

  double                 penetrance_indiv = pen.get_polygenic_penetrance(indiv);

  SL2_4 *= penetrance_indiv;

  // Siblings

  double trans_polygeno_sib = 0.0;
  double trans_geno_sib     = 0.0;
  double penetrance_sib     = 0.0;

  log_double posterior_sib           (0.0);
  log_double genotype_loop_total     (0.0);    
  log_double polygenotype_loop_total (0.0);

  const member_type& member = *indiv.member;

  member_type::sibling_const_iterator sib_loop = member.sibling_begin();

  for( ; sib_loop != member.sibling_end(); ++sib_loop)
  {
    polygenotype_loop_total = 0.0;
    for(size_t polygenotype_loop = 0; polygenotype_loop < mod.fpmm_sub_model.max_pgt(); ++polygenotype_loop)
    {
      trans_polygeno_sib  = trans_poly_calc.prob(polygenotype_loop,
                            mother.polygenotype,father.polygenotype,
                            mod.fpmm_sub_model.loci());

      if(!trans_polygeno_sib)
        continue;

      genotype_loop_total = 0.0;
      for(genotype_index genotype_loop = index_AA; genotype_loop != index_INVALID; ++genotype_loop)
      {
        trans_geno_sib       = mod.transm_sub_model.prob(genotype_loop,mother.genotype,father.genotype);

        if(!trans_geno_sib)
          continue;

        sibling              = penetrance_info(*sib_loop,
                                   genotype_loop,polygenotype_loop);
        penetrance_sib       = pen.get_polygenic_penetrance(sibling);
        posterior_sib        = peeler_pointer->posterior(*sib_loop,
                               genetic_info(genotype_loop,polygenotype_loop));
        genotype_loop_total += trans_geno_sib * penetrance_sib * posterior_sib;
      }
      polygenotype_loop_total += trans_polygeno_sib * genotype_loop_total;
    }
    SL2_4 *= polygenotype_loop_total;
  }

  return SL2_4;
}

//===================================================================
//
//  FPMM_SL::int_nonfounder_SL6(...)
//
//===================================================================
log_double
FPMM_SL::int_nonfounder_SL6(const penetrance_info&          indiv,
		       	    const penetrance_info&          spouse) const
{
  // For definition of SL6, see equation #40.

  log_double SL6(1.0);

  penetrance_info        child;

  // Calculate psi and penetrance values for first portion of SL2:

  double trans_geno_child     = 0.0;
  double trans_polygeno_child = 0.0;
  double penetrance_child     = 0.0;

  log_double posterior_child         (0.0);
  log_double genotype_loop_total     (0.0);
  log_double polygenotype_loop_total (0.0);

  // Child loop begins here

  const member_type& member = *indiv.member;
  const member_type& mate   = *spouse.member;

  member_type::offspring_const_iterator child_loop = member.offspring_begin(mate);

  for( ; child_loop != member.offspring_end(); ++child_loop)
  {
    polygenotype_loop_total = 0.0;
    for(size_t polygenotype_loop = 0; polygenotype_loop < mod.fpmm_sub_model.max_pgt(); ++polygenotype_loop)
    {
      trans_polygeno_child = trans_poly_calc.prob(polygenotype_loop,
                             indiv.polygenotype,spouse.polygenotype,
                             mod.fpmm_sub_model.loci());

      if(!trans_polygeno_child)
        continue;

      genotype_loop_total = 0.0;
      for(genotype_index genotype_loop = index_AA; genotype_loop != index_INVALID; ++genotype_loop)
      {
        trans_geno_child     = mod.transm_sub_model.prob(genotype_loop,indiv.genotype,spouse.genotype);

        if(!trans_geno_child)
          continue;

        child                = penetrance_info(*child_loop,
                               genotype_loop,polygenotype_loop);
        penetrance_child     = pen.get_polygenic_penetrance(child);
        posterior_child      = peeler_pointer->posterior(*child_loop,
                               genetic_info(genotype_loop,polygenotype_loop));
        genotype_loop_total += trans_geno_child * penetrance_child * posterior_child;
      }
      polygenotype_loop_total += trans_polygeno_child * genotype_loop_total;
    }
    SL6 *= polygenotype_loop_total;
  }
  return SL6;
}

}}
