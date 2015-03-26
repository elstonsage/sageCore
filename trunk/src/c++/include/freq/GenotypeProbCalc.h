#ifndef FREQ_GENOTYPE_PROB_H
#define FREQ_GENOTYPE_PROB_H

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
#include "freq/AlleleRemapping.h"
#include "freq/RemappedAlleleFreqs.h"

namespace SAGE {
namespace FREQ {

class GenotypeProbCalc
{
  public:

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    explicit GenotypeProbCalc(const MLOCUS::genotype_model & gmodel)
    : my_gmodel(gmodel),
      my_use_inbreeding(false),
      my_last_inbreeding_value(QNAN)
    {
      my_cache.clear();
      
      my_cache.resize(my_gmodel.unphased_genotype_count());
    }

    ///
    /// Copy constructor.
    GenotypeProbCalc(const GenotypeProbCalc & other) : my_gmodel(other.my_gmodel), my_use_inbreeding(other.my_use_inbreeding), my_cache(other.my_cache) { }
    
  //@}
  
  /// @name Misc.
  //@{
    
    static double getInbredHomozygousFrequency(double freq, double inbreeding) 
    {
      return (freq * freq) + freq * (1.0 - freq) * inbreeding;
    }

    static double getHomozygousFrequency(double freq) 
    {
      return (freq * freq);
    }

    static double getInbredHeterozygousFrequency(double freq1, double freq2, double inbreeding) 
    {
      return 2.0 * freq1 * freq2 * (1.0 - inbreeding);
    }

    static double getHeterozygousFrequency(double freq1, double freq2) 
    {
      return 2.0 * freq1 * freq2;
    }
    
    ///
    /// Populates the cache with genotype probabilities.
    /// \param mgr The ParameterMgr containing, among other things, the inbreeding coefficient (if used).
    /// \param allele_freqs The remapped allele frequencies
    void populateCache(const MAXFUN::ParameterMgr & mgr, const RemappedAlleleFreqs & allele_freqs)
    {
      if(my_use_inbreeding)
      {
        double inbreeding = mgr.getParameter("Inbreeding coeff.", "Inbreeding coeff.").getCurrentEstimate();

        for(MLOCUS::unphased_genotype_iterator i = my_gmodel.unphased_genotype_begin();
            i != my_gmodel.unphased_genotype_end(); ++i)
        {
          double f1 = allele_freqs.getFrequency(i->allele1().id()),
                 f2 = allele_freqs.getFrequency(i->allele2().id());
                 
          my_cache[i->get_id()] = i->homozygous() ? getInbredHomozygousFrequency(f1, inbreeding) :
                                                    getInbredHeterozygousFrequency(f1, f2, inbreeding);
        }
        my_last_inbreeding_value = inbreeding;
      }  
      else // No inbreeding
      { 
        for(MLOCUS::unphased_genotype_iterator i = my_gmodel.unphased_genotype_begin();
            i != my_gmodel.unphased_genotype_end(); ++i)
        {
          double f1 = allele_freqs.getFrequency(i->allele1().id()),
                 f2 = allele_freqs.getFrequency(i->allele2().id());

          my_cache[i->get_id()] = i->homozygous() ? getHomozygousFrequency(f1) :
                                                    getHeterozygousFrequency(f1, f2);
        }
      }  
    }

    ///
    /// Populates the cache with genotype probabilities.
    ///
    /// Note: This version of populateCache extracts allele frequencies directly from the ParameterMgr,
    ///       and does not take allele remapping into account.
    ///
    /// \param mgr The ParameterMgr containing the allele frequencies
    void populateCache(const MAXFUN::ParameterMgr & mgr)
    {
      if(my_use_inbreeding)
      {
        double inbreeding = mgr.getParameter("Inbreeding coeff.", "Inbreeding coeff.").getCurrentEstimate();

        for(MLOCUS::unphased_genotype_iterator i = my_gmodel.unphased_genotype_begin();
            i != my_gmodel.unphased_genotype_end(); ++i)
        {
          double f1 = mgr.getParameter(i->allele1().id()).getCurrentEstimate(),
                 f2 = mgr.getParameter(i->allele2().id()).getCurrentEstimate();

          my_cache[i->get_id()] = i->homozygous() ? getInbredHomozygousFrequency(f1, inbreeding) :
                                                    getInbredHeterozygousFrequency(f1, f2, inbreeding);
        }
        my_last_inbreeding_value = inbreeding;
      }  
      else // No inbreeding
      {   
        for(MLOCUS::unphased_genotype_iterator i = my_gmodel.unphased_genotype_begin();
            i != my_gmodel.unphased_genotype_end(); ++i)
        {
          double f1 = mgr.getParameter(i->allele1().id()).getCurrentEstimate(),
                 f2 = mgr.getParameter(i->allele2().id()).getCurrentEstimate();

          my_cache[i->get_id()] = i->homozygous() ? getHomozygousFrequency(f1) :
                                                    getHeterozygousFrequency(f1, f2);
        }
      }  
    }

    ///
    /// Returns the genotype probability.
    /// \param g The genotype id
    double getProb(uint g) const { 
    //cout << g << " -> " << my_cache[g] << endl;
    return my_cache[g]; }
    
    ///
    /// Sets whether or not to use inbreeding.
    void setUseInbreeding(bool i) { my_use_inbreeding = i; }
    
    ///
    /// Returns the genotype model bound to this genotype prob calc.
    const MLOCUS::genotype_model & getGmodel() const { return my_gmodel; }
    
    /// Returns \c true if inbreeding is in use and the inbreeding value has changed
    /// since the last maximization iteration. \c false otherwise.
    ///
    /// \param mgr The ParameterMgr being maximized.
    bool has_inbreeding_changed(const MAXFUN::ParameterMgr& mgr)
    {
      if(!my_use_inbreeding) return false;
      
      double inbreeding = mgr.getParameter("Inbreeding coeff.", "Inbreeding coeff.").getCurrentEstimate();

      return my_last_inbreeding_value != inbreeding;
    }

    ///
    /// Debugging only!
    void dump() const
    {
      std::cout << "GenotypeProbCalc dump:" << std::endl;
      
      for(MLOCUS::unphased_genotype_iterator i = my_gmodel.unphased_genotype_begin();
          i != my_gmodel.unphased_genotype_end(); ++i)
      {
        std::cout << "Genotype " << i->get_id() << " (" << i->name() << ") " << my_cache[i->get_id()] << std::endl;
      }
      
      std::cout << std::endl;
    }

  //@}
  
  private:

    GenotypeProbCalc& operator= (const GenotypeProbCalc & other); // Disallowed
    
    const MLOCUS::genotype_model & my_gmodel;

    bool my_use_inbreeding;
    
    double my_last_inbreeding_value;

    std::vector<double> my_cache;
};
  
} // End namespace FREQ
} // End namespace SAGE

#endif
