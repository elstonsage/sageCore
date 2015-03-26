#include "tdtex/TransmissionTables.h"

namespace SAGE  {
namespace TDTEX {

//================================================
//
// AlleleTransmissionTable CONSTRUCTOR
//
//================================================
AlleleTransmissionTable::AlleleTransmissionTable(const MLOCUS::genotype_model& gmodel) :
  TransmissionTable (gmodel.allele_count()),
  my_gmodel         (gmodel)
{ }

//================================================
//  
// AlleleTransmissionTable COPY CONSTRUCTOR   
//
//================================================
AlleleTransmissionTable::AlleleTransmissionTable(const AlleleTransmissionTable & other) :
  TransmissionTable (other),
  my_gmodel         (other.my_gmodel)
{ }
                
//================================================
//  
// AlleleTransmissionTable operator=
//
//================================================
AlleleTransmissionTable& 
AlleleTransmissionTable::operator=(const AlleleTransmissionTable & other)
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
//  clone()
//
//================================================
AlleleTransmissionTable *
AlleleTransmissionTable::clone() const
{
  return new AlleleTransmissionTable(*this);
}

//================================================
//
//  informative(...)
//
//================================================
bool
AlleleTransmissionTable::informative(const TransmissionList& xmit)
{
  return xmit.get_list().size();
}

//================================================
//
//  count_transmission(...)
//
//================================================
void
AlleleTransmissionTable::count_transmissions(const TransmissionList& xmit)
{
  // Loop across all the transmissions in the list and increment each count:
  for(TransmissionVector::const_iterator x = xmit.get_list().begin(); x != xmit.get_list().end(); ++x)
  {
    get_counts()(x->transmitted.id(), x->not_transmitted.id()) += 1;
  }
}

//================================================
//
//
//
//================================================
string
AlleleTransmissionTable::heading(size_t n) const
{
  return (my_gmodel.allele_begin() + n)->name();
}

//================================================
//
//
//
//================================================
string
AlleleTransmissionTable::units() const
{
  return "Allele";
}


} // End namespace TDTEX
} // End namespace SAGE
