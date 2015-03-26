#include "tdtex/Transmission.h"

namespace SAGE  {
namespace TDTEX {

//===============================================================
//
// transmission_intersection(...)
//
//===============================================================
TransmissionList 
transmission_intersection(const TransmissionList& t1, const TransmissionList& t2_orig)
{
  // Create TransmissionList:
  TransmissionList t;

  // If there are errors in either constituent TransmissionList, then label
  // this TransmissionList as in error:
  if(t1.get_error() || t2_orig.get_error())
  {
    t.set_error(true);
    return t;
  }

  // Copy over 2nd-person transmission list:
  TransmissionList t2(t2_orig);

  // For each Transmission entry in the t1 list, find the corresponding
  // Transmission in the t2 list:

  for(TransmissionVector::const_iterator i = t1.get_list().begin(); i != t1.get_list().end(); ++i)
  {
    // Find matching transmission from t2
    TransmissionVector::iterator j = std::find(t2.get_list().begin(), t2.get_list().end(), *i);

    // If no match, keep looking
    if(j == t2.get_list().end())
      continue;

    // Otherwise, copy the transmission, remove the matching one from t2
    Transmission strans = *i;
    t2.get_list().erase(j);

    // Forget the source, if the matching transmission does not have one
    if(!j->source)
      strans.source = NULL;

    // Store the new transmission
    t.get_list().push_back(strans);
  }

  return t;
}

void
TransmissionList::dump() const
{
  std::cout << "Transmissionlist: " << *this << std::endl;
}

//===============================================================
//
// Output: TransmissionList
//
//===============================================================
std::ostream& operator<<(std::ostream& o, const TransmissionList& t)
{
  o << '[';

  for(TransmissionVector::const_iterator i = t.get_list().begin(); i != t.get_list().end(); ++i)
  {
    if(i != t.get_list().begin())
      o << ", ";

    o << '(' << i->transmitted.name() << ',' << i->not_transmitted.name() << ')' << " from ";

    if(i->source && i->source->name().size())
      o << "parent '" << i->source->name() << '\'';
    else
      o << "unknown parent";
  }

  o << ']';

  return o;
}


} // End namespace TDTEX
} // End namespace SAGE
