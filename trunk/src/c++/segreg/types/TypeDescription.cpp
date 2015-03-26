#include "segreg/types/TypeDescription.h"

namespace SAGE   {
namespace SEGREG {

OUTPUT::Section convert_to_output(const TypeDescription& t)
{
  OUTPUT::Section s(t.get_name() + " : " + t.get_description());
  
  OUTPUT::Table basics;
  
  basics << OUTPUT::TableColumn("State ID")
         << OUTPUT::TableColumn("State Description");

  for(TypeDescription::StateIterator i = t.begin(); i != t.end(); ++i)
  {
    basics << (OUTPUT::TableRow() << i->get_index() << i->get_name());
  }

  s << basics;
  
  return s;
}

} // End Namespace SEGREG
} // End Namespace SAGE


