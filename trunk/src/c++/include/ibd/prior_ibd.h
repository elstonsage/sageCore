#ifndef __PRIOR_IBD_H
#define __PRIOR_IBD_H

//=======================================================================
// File:   prior_ibd.cpp
//
// Notes:  Compute the prior probability of allele sharing IBD
//         given relationship
//
// Copyright(c) 1998 RC Elston
//  All Rights Reserved
//=======================================================================

#include "ibd/definitions.h"

namespace SAGE {

class RefPriorIBD
{
  public:

    RefPriorIBD();
    RefPriorIBD(const SAGE::RPED::RefPedigree& rp);
    RefPriorIBD(const SAGE::FPED::Pedigree&    fp);

    void compute(const SAGE::RPED::RefPedigree& rp);
    void compute(const SAGE::FPED::Pedigree&    fp);

    void clear()
    {
      prior_ibd.clear();
      size = 0;
    }
      
    // P(i1 & i2 share n alleles IBD)
    double f(size_t i1, size_t i2, size_t n) const
    {
      size_t x1 = std::min(i1,i2);
      size_t x2 = std::max(i1,i2);

      if( !prior_ibd.size() || x2 >= size || n > 2 )
        return std::numeric_limits<double>::quiet_NaN();
      
      double f0, f2;
      
      if( x1 == x2 ) 
      {
        f0 = 0.0;
        f2 = 1.0;
      }       
      else
      {
        f0 = prior_ibd[ x2*size + x1 ];
        f2 = prior_ibd[ x1*size + x2 ];
      }
      switch(n)
      {
        case  0 : return f0;
        case  1 : return 1.0 - f0 - f2;
        case  2 : return f2;
        default : return std::numeric_limits<double>::quiet_NaN();
      }
    }

    // P(i1 & i2 share n alleles IBD)
    void f(size_t i1, size_t i2, double &f0, double &f1, double &f2) const
    {
      size_t x1 = std::min(i1,i2);
      size_t x2 = std::max(i1,i2);

      if( !prior_ibd.size() || x2 >= size )
      {
        f0 = f1 = f2 = std::numeric_limits<double>::quiet_NaN();
        return;
      }
      
      if( x1 == x2 ) 
      {
        f0 = 0.0;
        f1 = 0.0;
        f2 = 1.0;
      }       
      else
      {
        f0 = prior_ibd[ x2*size + x1 ];
        f2 = prior_ibd[ x1*size + x2 ];
        f1 = 1.0 - f0 - f2;
      }
    }

  protected:

    typedef vector<double> prior_matrix;
    
    prior_matrix    prior_ibd;
    size_t          size;
};

}

#endif
