#ifndef LIKELIHOOD_ELEMENTS_H
#define LIKELIHOOD_ELEMENTS_H

#include "segreg/PedigreeDataSet.h"
#include "segreg/model.h"
#include "segreg/RegPenetranceCalculator.h"
#include "segreg/types/TypeDescription.h"
#include "boost/bind.hpp"

namespace SAGE {
namespace SEGREG {
  
class LikelihoodElements
{
  public:
  
    LikelihoodElements(const FPED::Multipedigree& ped_data,
                       const model&               mdl,
                       bool                       use_ascertainment);
                       
    int update();
    
    const TypeDescription& get_type_description() const;
    
    double get_frequency    (const TypeDescription::State&         state) const;
    double get_frequency    (const TypeDescription::StateIterator& state) const;
    double get_transmission (const TypeDescription::State&         istate,
                             const TypeDescription::State&         mstate,
                             const TypeDescription::State&         fstate)    const;
    double get_transmission (const TypeDescription::StateIterator& istate,
                             const TypeDescription::StateIterator& mstate,
                             const TypeDescription::StateIterator& fstate)    const;
    
    PenetranceContext get_penetrance_context() const;
    
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::State&         state) const;
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::StateIterator& state) const;

    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::State&         state,
                          const PenetranceContext&              context)  const;
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::StateIterator& state,
                          const PenetranceContext&              context)  const;
  
    const RegPenetranceCalculator& get_pc() const { return my_pen_calc; } 
  
  private:
  
    void build_type_description();
  
    boost::function<double (const genotype_index&)>  my_freq_function;
    boost::function<double (const genotype_index&,
                            const genotype_index&,
                            const genotype_index&)>  my_transm_function;
                            
    RegPenetranceCalculator my_pen_calc;
    TypeDescription         my_type_description;
};

}
}

#include "segreg/LikelihoodElements.ipp"

#endif

