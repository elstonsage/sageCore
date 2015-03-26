#ifndef LODLINK_TRANS_CALCULATOR_H
#define LODLINK_TRANS_CALCULATOR_H
//============================================================================
// File:      trans_calculator.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/14/2 - created.                                   djb
//                                                                          
// Notes:     defines a class to calculate transition probabilities.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <ostream>
#include "lodlink/mle_sub_model.h"
#include "lodlink/definitions.h"

namespace SAGE
{

namespace LODLINK
{

class trans_calculator
{
  public:
    trans_calculator(const mle_sub_model& mle);  
  
    const mle_sub_model&  mle() const;
    
    double transition(const joint_genotype& mom, 
                      const joint_genotype& dad,
                      const joint_genotype& kid ) const;
                      
  private:
    static bool  kd(const MLOCUS::allele& one, const MLOCUS::allele& two);    // Kronecker delta.
    double  transmission(const joint_genotype& jg, LODLINK::sex s, const haplotype& h) const;
    
    trans_calculator(const trans_calculator& other);
    trans_calculator&  operator=(const trans_calculator& other);
  
    //Data members.
    const mle_sub_model&  my_mle;
};

#include "lodlink/trans_calculator.ipp"
}
}

#endif


