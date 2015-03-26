#ifndef FREQ_TRANSMISSION_PROB_H
#define FREQ_TRANSMISSION_PROB_H

#include <vector>
#include <list>
#include "error/errormanip.h"
#include "error/errorstream.h"
#include "numerics/log_double.h"
#include "peeling/peeler3.h"
#include "util/AutoTrace.h"
#include "rped/loop.h"   
#include "rped/rped.h"
#include "maxfunapi/maxfunapi.h"
#include "freq/Sample.h"
#include "freq/Cache.h"

namespace SAGE {
namespace FREQ {

class TransmissionProbCalc
{
  public:

  /// @name Constructors
  //@{
  
    TransmissionProbCalc() { my_cache.clear(); }
    
    TransmissionProbCalc(const TransmissionProbCalc & other) : my_cache(other.my_cache) { }

    TransmissionProbCalc& operator=(const TransmissionProbCalc & other) { if(this != &other) { my_cache = other.my_cache; } return *this; }
    
  //@}
  
  /// @name Misc.
  //@{
    
    ///
    /// Populates the cache with genotype transmission probabilities.
    void populateCache(const MLOCUS::genotype_model & gmodel)
    {
      my_cache.clear();
      
      size_t genotype_count = gmodel.unphased_genotype_count();
      
      my_cache.resize(genotype_count, std::vector<std::vector<double> > (genotype_count, std::vector<double> (genotype_count, 0.0)));
      
      for(MLOCUS::unphased_genotype_iterator g_p1 = gmodel.unphased_genotype_begin();
          g_p1 != gmodel.unphased_genotype_end(); ++g_p1)
      {
        for(MLOCUS::unphased_genotype_iterator g_p2 = gmodel.unphased_genotype_begin();
            g_p2 != gmodel.unphased_genotype_end(); ++g_p2)
        {
          MLOCUS::allele p1_a1 = g_p1->allele1(),
                         p1_a2 = g_p1->allele2(),
                         p2_a1 = g_p2->allele1(),
                         p2_a2 = g_p2->allele2();

          my_cache[g_p1->get_id()][g_p2->get_id()][gmodel.get_unphased_genotype(p1_a1, p2_a1).get_id()] += 0.25;
          my_cache[g_p1->get_id()][g_p2->get_id()][gmodel.get_unphased_genotype(p1_a1, p2_a2).get_id()] += 0.25;
          my_cache[g_p1->get_id()][g_p2->get_id()][gmodel.get_unphased_genotype(p1_a2, p2_a1).get_id()] += 0.25;
          my_cache[g_p1->get_id()][g_p2->get_id()][gmodel.get_unphased_genotype(p1_a2, p2_a2).get_id()] += 0.25;
        }
      }
    }

    ///
    /// Returns the transmissions probability.
    /// \param gchild The child's genotype id
    /// \param gparent1 The parent 1's genotype id
    /// \param gparent2 The parent 2's genotype id
    double getProb(uint gchild, uint gparent1, uint gparent2) const
    {
      return my_cache[gparent1][gparent2][gchild];
    }

  //@}
  
  private:

    /// Index #1 is parent1
    /// Index #2 is parent2
    /// Index #3 is child
    typedef std::vector<std::vector<std::vector<double> > > ProbCache;
    
    ProbCache my_cache;
};
  
} // End namespace FREQ
} // End namespace SAGE

#endif
