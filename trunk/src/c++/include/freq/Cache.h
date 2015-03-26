#ifndef FREQ_CACHE
#define FREQ_CACHE

#include <vector>
#include <map>
#include "numerics/log_double.h"
#include "peeling/cache3.h"
#include "fped/fped.h"

namespace SAGE    {
namespace peeling {

template<>
class individual_cache<MLOCUS::unphased_genotype, log_double>
{
  public:

  /// @name Typedefs
  //@{
  
    typedef MLOCUS::unphased_genotype data_type;
    typedef log_double                result_type;
    
    typedef FPED::FilteredMultipedigree::member_type member_type;
    
    typedef std::vector<bool> CachedVector;
    
    typedef std::vector<log_double> ValVector;
    
    struct MateLookup
    {
      MateLookup() { }

      MateLookup(size_t _mate_id, size_t _idx) : mate_id(_mate_id), idx(_idx) { }

      MateLookup(const MateLookup & other) : mate_id(other.mate_id), idx(other.idx) { }

      MateLookup& operator= (const MateLookup & other) { if(this != & other) { mate_id = other.mate_id; idx = other.idx; } return *this; }
      
      size_t mate_id; // The id of the mate
      size_t idx;     // The beginning idx in the vector
    };
    
    typedef std::vector<MateLookup> MateLookupVector;

  //@}
  
  /// @name Constructors / operators
  //@{

    ///
    /// Constructor.
    individual_cache() { }
    
    ///
    /// Copy constructor.
    individual_cache(const individual_cache & other) : my_cached(other.my_cached), my_vals(other.my_vals), my_mate_lookups(other.my_mate_lookups) { }

    ///
    /// Assignment operator.
    individual_cache& operator=(const individual_cache & other) 
    { 
      if(this != &other) 
      { 
        my_cached       = other.my_cached;
        my_vals         = other.my_vals; 
        my_mate_lookups = other.my_mate_lookups;
      } 
      return *this; 
    }
    
  //@}
  
  /// @name Misc
  //@{
  
    ///
    /// Sets up the necessary vectors / structs for caching values.
    /// \param genotype_count Number of genotypes
    /// \param ind The individual this cache represents
    void setup(size_t genotype_count, const member_type & ind)
    {
      my_cached       . resize (genotype_count, false);
      my_vals         . resize (genotype_count);
      my_mate_lookups . clear  ();
      
      for(FPED::MateConstIterator mate_itr = ind.mate_begin(); mate_itr != ind.mate_end(); ++mate_itr)
      {
        my_mate_lookups.push_back(MateLookup(mate_itr->mate().mpindex(), my_vals.size()));

        my_vals   .resize(my_vals   .size() + genotype_count);
        my_cached .resize(my_cached .size() + genotype_count, false);
      }
      
      clearCachedVals();
    }
  
    ///
    /// Clears any cached vals, but maintains the structure of the cache.
    void clearCachedVals()
    {
      my_cached.assign(my_cached.size(), false);
    }
  
  //@}

  /// @name Anterior
  //@{
    
    ///
    /// Anterior cached?
    bool anterior_cached(const data_type & g) const { return my_cached[g.get_id()]; }
    
    ///
    /// Anterior.
    const result_type & anterior(const data_type & g) const { return my_vals[g.get_id()]; }
    
    ///
    /// Non-const anterior
    result_type & anterior(const data_type & g)
    {
      my_cached[g.get_id()] = true;

      return my_vals[g.get_id()];
    }
    
  //@}
  
  /// @name Posterior
  //@{

    ///
    /// Posterior-with-mate cached?
    bool posterior_with_mate_cached (const member_type & mate, const data_type & g) const
    {
      for(MateLookupVector::const_iterator i = my_mate_lookups.begin(); i != my_mate_lookups.end(); ++i)
      {
        if(i->mate_id == mate.mpindex())
        {
          return my_cached[i->idx + g.get_id()];
        }
      }
      
      throw std::exception();
    }

    ///
    /// Posterior with mate.
    const result_type & posterior_with_mate (const member_type & mate, const data_type &g) const
    {
      for(MateLookupVector::const_iterator i = my_mate_lookups.begin(); i != my_mate_lookups.end(); ++i)
      {
        if(i->mate_id == mate.mpindex())
        {
          return my_vals[i->idx + g.get_id()];
        }
      }
      
      throw std::exception();
    }

    ///
    /// Non-const posterior-with-mate
    result_type & posterior_with_mate (const member_type & mate, const data_type & g)
    {
      for(MateLookupVector::iterator i = my_mate_lookups.begin(); i != my_mate_lookups.end(); ++i)
      {
        if(i->mate_id == mate.mpindex())
        {
          size_t real_idx = i->idx + g.get_id();
          
          my_cached[real_idx] = true;

          return my_vals[real_idx];
        }
      }
      
      throw std::exception();
    }

    
  //@}
  
  /// @name Debug
  //@{
  
    void dump() const;
    
  //@}

  private:

  /// An explanation of the CachedValVector:
  ///
  /// The vector is a single-dimensional storage of a bunch of different pieces of cached information.
  //
  /// Indices 0 through genotype_count - 1 : The anterior value for the given index
  ///
  /// Indices genotype_count + (mate_num * genotype_count) through
  /// genotype_count + (mate_num * genotype_count) + genotype_count - 1: The posterior_with_mate for the given mate and genotype.
  
    CachedVector     my_cached;
    ValVector        my_vals;
    MateLookupVector my_mate_lookups;

};

} // End namespace FREQ
} // End namespace SAGE

#endif
