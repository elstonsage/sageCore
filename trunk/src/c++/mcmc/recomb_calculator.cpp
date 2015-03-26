//==========================================================================
//  File:     recomb_calculator.h                                           
//                                                                          
//  Author:   Geoff Wedig                                                   
//                                                                          
//  History:  0.1 Initial Implementation                                    
//            1.0 Updated to new libraries                       yjs May. 04
//                                                                          
//  Notes:    Calculates the inheritance probabilities of intervals         
//            between markers based upon Lander-Green Inheritance patterns  
//            observed at those markers.  It can also claculate the change  
//            in ratio when those patterns change.                          
//                                                                          
//  Copyright (c) 2004 R. C. Elston                                         
//  All Rights Reserved                                                     
//==========================================================================

#include "mcmc/recomb_calculator.h"

namespace SAGE
{

namespace MCMC
{

recombination_calculator::recombination_calculator(const pedigree_region&    r,
                                                   const mcmc_data_accessor& d)
                        : my_ped_region(r), my_data(d)
{
  my_log_thetas.resize(r.get_region().locus_count()-1);
  my_log_one_minus_thetas.resize(r.get_region().locus_count()-1);
  my_theta_ratios.resize(r.get_region().locus_count()-1);

  for( size_t i = 0; i < r.get_region().locus_count()-1; ++i )
  {
    double theta = r.get_region().locus(i).locus_theta(1);

    // This is BAD, but needed right now.
    if( theta <= 0.0 ) theta = numeric_limits<double>::epsilon();
    if( theta >= 1.0 ) theta = (double) 1.0 - numeric_limits<double>::epsilon();

    const double lth  = log(theta);

#ifdef WINDOWS
    // g++ doesn't have log1p *sigh*
    const double lmth = log(1.0-theta);
#else
    const double lmth = log1p(-theta);
#endif

    my_log_thetas[i]           = lth;
    my_log_one_minus_thetas[i] = lmth;
    my_theta_ratios[i]         = lth - lmth;
  }
}

void
recombination_calculator::dump_recombination_calculator(ostream& o) const
{
  o << endl << "recombination_calculator dump :" << endl;

  o << endl << "  log_thetas :" << endl;
  for( size_t i = 0; i < my_log_thetas.size(); ++i )
  {
    if( i && i % 5 == 0 )
      o << endl;

    o << "\t" << my_log_thetas[i];
  }
  o << endl;

  o << endl << "  log_one_minus_thetas :" << endl;
  for( size_t i = 0; i < my_log_one_minus_thetas.size(); ++i )
  {
    if( i && i % 5 == 0 )
      o << endl;

    o << "\t" << my_log_one_minus_thetas[i];
  }
  o << endl;

  o << endl << "  theta_ratios :" << endl;
  for( size_t i = 0; i < my_theta_ratios.size(); ++i )
  {
    if( i && i % 5 == 0 )
      o << endl;

    o << "\t" << my_theta_ratios[i];
  }
  o << endl;
}

} // end of namespace MCMC

} // end of namespace SAGE
