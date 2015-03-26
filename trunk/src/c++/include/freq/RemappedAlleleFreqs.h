#ifndef FREQ_REMAPPED_ALLELE_FREQS
#define FREQ_REMAPPED_ALLELE_FREQS

#include "maxfunapi/maxfunapi.h"

namespace SAGE {
namespace FREQ {

class RemappedAlleleFreqs
{
  public:
  
  /// @name Constructors
  //@{

    ///
    /// Constructor.
    explicit RemappedAlleleFreqs(size_t allele_count)
    {
      my_freqs.resize(allele_count);
    }
    
    ///
    /// Copy constructor.
    RemappedAlleleFreqs(const RemappedAlleleFreqs & other) : my_freqs(other.my_freqs) { }
    
    ///
    /// Operator=
    RemappedAlleleFreqs& operator= (const RemappedAlleleFreqs & other)
    {
      if(this != &other)
      {
        my_freqs = other.my_freqs;
      }
      
      return *this;
    }
    
  //@}

  /// @name Misc
  //@{
  
    /// Figures out the remapped allele frequencies.
    /// \returns True if any of the remapped alleles frequencies have changed, false otherwise
    bool update(const MAXFUN::ParameterMgr & mgr, const AlleleRemapping & allele_remapping)
    {
      bool remapped_freqs_changed = false;
      
      for(size_t remapped_id = 0; remapped_id < allele_remapping.getRemappedAlleleCount(); ++remapped_id)
      {
        size_t unmapped_id = allele_remapping.getUnmappedId(remapped_id);
        
        if(unmapped_id == (size_t)-1)
        {
          double remap_freq = 0.0;
          
          for(std::set<size_t>::const_iterator unmapped_id  = allele_remapping.getRemapSet().begin (); 
                                               unmapped_id != allele_remapping.getRemapSet().end   (); ++unmapped_id)
          {
            remap_freq += mgr.getParameter(*unmapped_id).getCurrentEstimate();
          }
          
          remapped_freqs_changed |= my_freqs[remapped_id] != remap_freq;

          my_freqs[remapped_id] = remap_freq;
        }
        else // Not a remapped allele
        {
          const MAXFUN::Parameter & param = mgr.getParameter(unmapped_id);
          
          my_freqs[remapped_id] = param.getCurrentEstimate();
          
          remapped_freqs_changed |= param.estimateChanged();
        }
      }
      
      return remapped_freqs_changed;
    }
    
    ///
    /// Returns the frequency of the given allele.
    /// \param allele_id The remapped allele id
    double getFrequency(size_t allele_id) const
    {
      return my_freqs[allele_id];
    }
    
    ///
    /// Debugging only!
    void dump() const
    {
      std::cout << "RemappedAlleleFreqs dump:" << std::endl;
      
      for(size_t i = 0; i < my_freqs.size(); ++i)
      {
        std::cout << "Remapped allele " << i << " = " << my_freqs[i] << std::endl;
      }
      
      std::cout << std::endl;
    }
    
  //@}
  
  
  private:
  
    std::vector<double> my_freqs;
};



} // End namespace FREQ
} // End namespace SAGE

#endif

