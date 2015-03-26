#ifndef AO_VALIDATOR_H
#define AO_VALIDATOR_H

#include "sampling/sampling.h"

namespace SAGE {
namespace AO   {

class Validator : public SAMPLING::IndividualValidator
{
  public:
    virtual bool isValid (size_t i, const SAMPLING::IndividualTraitData & trait_data) const;

    void setMultiPedigree(const FPED::Multipedigree& mp) { my_mp = &mp; }

  private:
    const FPED::Multipedigree* my_mp;
};

} // End namespace AO
} // End namespace SAGE

#endif
