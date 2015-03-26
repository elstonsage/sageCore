#include "tdtex/TransmissionTables.h"

namespace SAGE  {
namespace TDTEX {

//================================================
//
// CONSTRUCTOR
//
//================================================
GenotypeTransmissionTable::GenotypeTransmissionTable(const MLOCUS::genotype_model& gmodel) :
  TransmissionTable (gmodel.unphased_genotype_count()),
  my_gmodel         (gmodel)
{ }

//================================================
//
// COPY CONSTRUCTOR
//
//================================================
GenotypeTransmissionTable::GenotypeTransmissionTable(const GenotypeTransmissionTable & other) :
  TransmissionTable (other),
  my_gmodel         (other.my_gmodel)
{ }

//================================================
//
// operator=
//
//================================================
GenotypeTransmissionTable& 
GenotypeTransmissionTable::operator=(const GenotypeTransmissionTable & other)
{
  TransmissionTable::operator=(other);
  
  if(this != &other)
  {
    my_gmodel = other.my_gmodel;
  }
  
  return *this;
}

//================================================
//
// clone()
//
//================================================
GenotypeTransmissionTable * 
GenotypeTransmissionTable::clone() const
{
  return new GenotypeTransmissionTable(*this);
}


//================================================
//
//  informative()
//
//================================================
bool
GenotypeTransmissionTable::informative(const TransmissionList& xmit)
{
  return xmit.get_list().size() == 2;
}

//================================================
//
//  count_transmission(...)
//
//================================================
void
GenotypeTransmissionTable::count_transmissions(const TransmissionList& xmit)
{
  // Make sure this is a valid genotype transmission list (2 transmissions):
  if(xmit.get_list().size() != 2)
    return;

  // Fetch the allele ids:
  MLOCUS::allele trans_allele1     = xmit.get_list()[0].transmitted,
                 not_trans_allele1 = xmit.get_list()[0].not_transmitted,
                 trans_allele2     = xmit.get_list()[1].transmitted,
                 not_trans_allele2 = xmit.get_list()[1].not_transmitted;

  // Increment the proper allele count:
  get_counts()(my_gmodel.get_unphased_genotype(trans_allele1,     trans_allele2)     .get_id(),
               my_gmodel.get_unphased_genotype(not_trans_allele1, not_trans_allele2) .get_id()) += 1;
}

//================================================
//
//
//
//================================================
string
GenotypeTransmissionTable::heading(size_t n) const
{
  return my_gmodel.get_unphased_genotype(n).name();
}

//================================================
//
//
//
//================================================
string
GenotypeTransmissionTable::units() const
{
  return "Genotype";
}


} // End namespace TDTEX
} // End namespace SAGE
