//============================================================================
// File:      peeler.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/2/2 created        -djb
//                                                                          
// Notes:     inline implementation of peeler class.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  peeler
//============================================================================
//
inline const peeler::subped_type&
peeler::subpedigree() const
{
  return my_subpedigree;
}

inline size_t
peeler::trait() const
{
  return my_trait;
}

inline size_t
peeler::marker() const
{
  return my_marker;
}

inline const trans_calculator&
peeler::tcalc() const
{
  return my_tcalc;
}

// - Sum likelihood of and individual and his posterior over his phenoset.
//
inline log_double
peeler::sum_ind_and_posterior(const joint_genotype& mjg, const joint_genotype& fjg,
                              const member_type& ind)
{
  log_double  r(0);
  phenoset    ph_set(my_subpedigree, my_trait, my_marker, ind);
  
  for(phenoset_iterator iter = ph_set.begin(); iter != ph_set.end(); ++iter)
  {
    joint_pen_iter  jpi = *iter;
    double  trans = my_tcalc.transition(mjg, fjg, joint_genotype(jpi));
    if(trans == 0)
    {
      continue;
    }
    
      r += jpi.penetrance() * posterior(ind, jpi) * trans;
  }
  
  return r;
}



