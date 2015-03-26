//============================================================================
// File:      trans_calculator.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/14/2 - created.                                djb
//                                                                          
// Notes:     inlines for trans_calculator class.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  trans_calculator
//============================================================================
//
inline  
trans_calculator::trans_calculator(const mle_sub_model& mle)
    : my_mle(mle)
{}

inline const mle_sub_model&
trans_calculator::mle() const
{
  return my_mle;
}

// - Kronecker delta.
//
inline bool
trans_calculator::kd(const MLOCUS::allele& one, const MLOCUS::allele& two)
{
  return (one == two) ? true : false;
}

// - Calculate the probability that an individual of the given sex w. the given 
//   joint genotype transmits the given haplotype.  This equation is given in
//   the LODLINK 3.1 user documentation.
//
inline double
trans_calculator::transmission(const joint_genotype& jg, LODLINK::sex s, const haplotype& h) const
{
  MLOCUS::allele  t_allele1 = jg.tg.allele1();
  MLOCUS::allele  t_allele2 = jg.tg.allele2();
  MLOCUS::allele  m_allele1 = jg.mg.allele1();
  MLOCUS::allele  m_allele2 = jg.mg.allele2();
  
  if((t_allele1 != h.ta && t_allele2 != h.ta) ||
     (m_allele1 != h.ma && m_allele2 != h.ma)   )
  {
    return 0;
  }
  
  double  theta = my_mle.theta(s);
  
  double  term1 = kd(t_allele1, h.ta) && kd(m_allele1, h.ma);
  double  term2 = kd(t_allele2, h.ta) && kd(m_allele2, h.ma);
  double  term3 = kd(t_allele1, h.ta) && kd(m_allele2, h.ma);
  double  term4 = kd(t_allele2, h.ta) && kd(m_allele1, h.ma); 
  
  double  non_recomb = (1 - theta) * (term1 + term2);
  double  recomb     =      theta  * (term3 + term4);

  return  (non_recomb + recomb) / 2;
  
  /*
  // - This less readable version was an attempt to improve program speed.
  //   It ran at the same speed, however.
  //
  return  ((1 - theta) * (kd(t_allele1, h.ta) * kd(m_allele1, h.ma) + kd(t_allele2, h.ta) * kd(m_allele2, h.ma)) +
                theta  * (kd(t_allele1, h.ta) * kd(m_allele2, h.ma) + kd(t_allele2, h.ta) * kd(m_allele1, h.ma))) / 2;
  */              
}

// - Note: this differs from the LODLINK 3.1 user documentation, but per
//   GCW, it is necessary to prevent errors in the case of missing data because of the 
//   way phase is treated in mlocus.
//
inline double
trans_calculator::transition(const joint_genotype& mom,
                             const joint_genotype& dad,
                             const joint_genotype& kid ) const
{
  
  double  trans1 = transmission(mom, LODLINK::female, haplotype(kid.tg.allele1(), kid.mg.allele1()));
  double  trans2 = transmission(dad, LODLINK::male, haplotype(kid.tg.allele2(), kid.mg.allele2()));
  
  return trans1 * trans2;
  
  
  /* As described in 3.1 documentation.
  double  trans1 = transmission(mom, female, haplotype(kid.tg.allele1(), kid.mg.allele1()));
  double  trans2 = transmission(dad, male, haplotype(kid.tg.allele2(), kid.mg.allele2()));
  double  trans3 = transmission(mom, female, haplotype(kid.tg.allele2(), kid.mg.allele2()));
  double  trans4 = transmission(dad, male, haplotype(kid.tg.allele1(), kid.mg.allele1()));
  
  if(kid.tg.allele1() == kid.tg.allele2() &&
     kid.mg.allele1() == kid.mg.allele2()    )
  {
    return  trans1 * trans2;
  }
  else
  {
    return  trans1 * trans2 + trans3 * trans4;
  }
  */
}

