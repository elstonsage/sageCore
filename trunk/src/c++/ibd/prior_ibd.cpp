//=======================================================================
// prior_ibd.cpp
//
// Compute the prior probability of allele sharing IBD given relationship
//
// Copyright(c) 1998 RC Elston
//=======================================================================

#include <cmath>
#include "mped/mp.h"
#include "ibd/prior_ibd.h"

namespace SAGE
{

RefPriorIBD::RefPriorIBD()
{
  size = 0;
}

RefPriorIBD::RefPriorIBD(const SAGE::RPED::RefPedigree& rp)
{
  size = 0;
  compute(rp);
}

RefPriorIBD::RefPriorIBD(const SAGE::FPED::Pedigree& fp)
{
  size = 0;
  compute(fp);
}

void RefPriorIBD::compute(const SAGE::RPED::RefPedigree& rp)
{
  size = rp.member_count();
  prior_ibd.resize( size*size );

  RPED::RefPedigree::member_type const * mid1;
  RPED::RefPedigree::member_type const * mid2;

  int p1, p2, p3, p4;
  double p1s0, p1s1, p1s2, p2s0, p2s1, p2s2;
  for(size_t i = 0; i < size; ++i)
    for(size_t j = 0; j < i; ++j)
    {
      mid1=&rp.member_index(i);
      mid2=&rp.member_index(j);

      if(mid1->parent1() || mid2->parent1())
      {
        p1 = p2 = p3 = p4 = -1;

        if(mid1->parent1())
        {
          p1 = mid1->parent1()->index();
          p2 = mid1->parent2()->index();
        }

        if(mid2->parent1())
        {
          p3 = mid2->parent1()->index();
          p4 = mid2->parent2()->index();
        }

        // Make sure we compute the prob. in generational order
        if(max(max(p1, p2), max(p3, p4)) == max(p1, p2))
        {
          f(p1, j, p1s0, p1s1, p1s2);
          f(p2, j, p2s0, p2s1, p2s2);
        }
        else
        {
          f(p3, i, p1s0, p1s1, p1s2);
          f(p4, i, p2s0, p2s1, p2s2);
        }
      }
      else
      {
        p1s0 = p2s0 = 1;
        p1s1 = p2s1 = p1s2 = p2s2 = 0;
      }

#if DEBUG_PRIOR_IBD_COMPUTATION
      cout << mid1->name() << ' ' << mid2->name() << ' '
           << p1s0 << ' ' << p1s1 << ' ' << p1s2 << ' '
           << p2s0 << ' ' << p2s1 << ' ' << p2s2 << endl;
#endif

      // P(i,j,0)
      prior_ibd[i * size + j] = p1s0*p2s0 + p1s0*p2s1/2 + p1s1*p2s0/2
                              + p1s1*p2s1/4;
      // P(i,j,2)
      prior_ibd[j * size + i] = p1s2*p2s2 + p1s2*p2s1/2 + p1s1*p2s2/2
                              + p1s1*p2s1/4;
    }
}

void RefPriorIBD::compute(const SAGE::FPED::Pedigree& rp)
{
  size = rp.member_count();
  prior_ibd.resize( size*size );

  FPED::Member const * mid1;
  FPED::Member const * mid2;

  int p1, p2, p3, p4;
  double p1s0, p1s1, p1s2, p2s0, p2s1, p2s2;
  for(size_t i = 0; i < size; ++i)
    for(size_t j = 0; j < i; ++j)
    {
      mid1=&rp.member_index(i);
      mid2=&rp.member_index(j);

      if(mid1->parent1() || mid2->parent1())
      {
        p1 = p2 = p3 = p4 = -1;

        if(mid1->parent1())
        {
          p1 = mid1->parent1()->index();
          p2 = mid1->parent2()->index();
        }

        if(mid2->parent1())
        {
          p3 = mid2->parent1()->index();
          p4 = mid2->parent2()->index();
        }

        // Make sure we compute the prob. in generational order
        if(max(max(p1, p2), max(p3, p4)) == max(p1, p2))
        {
          f(p1, j, p1s0, p1s1, p1s2);
          f(p2, j, p2s0, p2s1, p2s2);
        }
        else
        {
          f(p3, i, p1s0, p1s1, p1s2);
          f(p4, i, p2s0, p2s1, p2s2);
        }
      }
      else
      {
        p1s0 = p2s0 = 1;
        p1s1 = p2s1 = p1s2 = p2s2 = 0;
      }

#if DEBUG_PRIOR_IBD_COMPUTATION
      cout << mid1->name() << ' ' << mid2->name() << ' '
           << p1s0 << ' ' << p1s1 << ' ' << p1s2 << ' '
           << p2s0 << ' ' << p2s1 << ' ' << p2s2 << endl;
#endif

      // P(i,j,0)
      prior_ibd[i * size + j] = p1s0*p2s0 + p1s0*p2s1/2 + p1s1*p2s0/2
                              + p1s1*p2s1/4;
      // P(i,j,2)
      prior_ibd[j * size + i] = p1s2*p2s2 + p1s2*p2s1/2 + p1s1*p2s2/2
                              + p1s1*p2s1/4;
    }
}

}
