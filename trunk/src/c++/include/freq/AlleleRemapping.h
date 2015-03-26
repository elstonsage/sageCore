#ifndef FREQ_ALLELE_REMAPPING
#define FREQ_ALLELE_REMAPPING

#include "maxfunapi/maxfunapi.h"

namespace SAGE {
namespace FREQ {

class AlleleRemapping
{
  public:
  
    ///
    /// Given the original (unremapped) and remapped gmodels, figures out the translation scheme between the two.
    void setup(const MLOCUS::genotype_model & original_gmodel, const MLOCUS::genotype_model & remapped_gmodel)
    {
      // Resize the vectors:
      my_ids.resize(remapped_gmodel.allele_count());

      // Populate the remap set:
      for(size_t i = 0; i < original_gmodel.allele_count(); ++i)
        my_remap_set.insert(i);
      
      // Loop through the alleles in the remapped model:
      for(MLOCUS::allele_iterator remapped_al = remapped_gmodel.allele_begin();
          remapped_al != remapped_gmodel.allele_end(); ++remapped_al)
      {
        std::string remapped_allele_name = remapped_al->name();
        
        size_t unmapped_id = 0;
        
        // If it's the remap allele, then set the remap id to -1 (which will be specially processed by the update function).
        if(remapped_allele_name == "~remap")
        {
          unmapped_id = (size_t)-1;
        }

        // If it's not a remap allele, then set the remap id and remove the remap id from the remap set.
        else
        {
          unmapped_id = original_gmodel.get_allele(remapped_allele_name).id();
          
          my_remap_set.erase(unmapped_id);
        }

        my_ids[remapped_al->id()] = unmapped_id;
      }
    }

    ///
    /// Returns the number of remapped alleles.    
    size_t getRemappedAlleleCount() const { return my_ids.size(); }

    ///
    /// Returns the unmapped id.
    size_t getUnmappedId(size_t remapped_id) const { return my_ids[remapped_id]; }

    ///
    /// Returns the set of remap alleles.
    const std::set<size_t> & getRemapSet() const { return my_remap_set; }

    void dump() const
    {
      std::cout << "AlleleRemapping dump:" << std::endl;
      
      for(size_t i = 0; i < my_ids.size(); ++i)
      {
        std::cout << "Remapped idx " << i << " maps to original id ";
        
        if(my_ids[i] == (size_t)-1)
        {
          for(std::set<size_t>::const_iterator j = my_remap_set.begin(); j != my_remap_set.end(); ++j)
            std::cout << *j << " ";
        }
        else
        {
          std::cout << my_ids[i];
        }
        
        std::cout << std::endl;
      }
      
      std::cout << std::endl;
    }

  private:
  
    std::vector<size_t> my_ids;
    
    std::set<size_t> my_remap_set;

};



} // End namespace FREQ
} // End namespace SAGE

#endif

