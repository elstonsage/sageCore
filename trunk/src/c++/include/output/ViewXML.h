#ifndef OUTPUT_VIEW_XML_H
#define OUTPUT_VIEW_XML_H

#include "output/Graph.h"
#include "output/Misc.h"
#include "output/RenderingRules.h"
#include "output/Table.h"
#include "util/OutlineCntr.h"

namespace SAGE   {
namespace OUTPUT {

class ViewXML
{
  public:

  /// @name Output formatting
  //@{
  
    ///
    /// Renders to a file.
    template<typename T>
    static void renderToFile(const T & element, const std::string & filename)
    {
      std::ofstream o;      
      o.open(filename.c_str());      

      o << render(element) << std::flush;      
    }

    template<typename T>
    static std::string render(const T & element)
    {
      // If parent_element is null, spawn a new document:
      try
      {
        xmlpp::Document document;

        xmlpp::Element * root_node = document.create_root_node("XML_DATA_FILE");

        internal_render(element, root_node);

        std::ostringstream s;
        
        s << document.write_to_string() << std::flush;
        
        return s.str();
      }
      catch(const std::exception& ex)
      {
        return "";
        // std::cout << "Exception caught: " << ex.what() << std::endl;
      }
    }
    
  private:
  
    template<typename T> 
    static void internal_render(const T & element, xmlpp::Element * parent_element)
    {
      // Create child element:
      xmlpp::Element * child_element = parent_element->add_child(element.getType());

      // Add attributes:
      typename T::AttributeMap attr_map;

      element.populateAttrMap(attr_map);
      
      for(typename T::AttributeMap::const_iterator i = attr_map.begin(); i != attr_map.end(); ++i)
        child_element->set_attribute(i->first, i->second);

      // Process sub-children:
      for(size_t i = 0; i < element.getVector().size(); ++i)
      {
        // Table, Double, Int, String, Graph, List, NamedDouble, NamedInt, NamedString, RenderingRules, Section, SpannedCell, TableColumn, TableRow, UnavailableCell
        
             if(element.getVector().template isType<Table>           (i)) internal_render(element.getVector().template getAbs<Table>           (i), child_element);
        else if(element.getVector().template isType<Double>          (i)) internal_render(element.getVector().template getAbs<Double>          (i), child_element);
        else if(element.getVector().template isType<Int>             (i)) internal_render(element.getVector().template getAbs<Int>             (i), child_element);
        else if(element.getVector().template isType<String>          (i)) internal_render(element.getVector().template getAbs<String>          (i), child_element);
        else if(element.getVector().template isType<Graph>           (i)) internal_render(element.getVector().template getAbs<Graph>           (i), child_element);
        else if(element.getVector().template isType<List>            (i)) internal_render(element.getVector().template getAbs<List>            (i), child_element);
        else if(element.getVector().template isType<NamedDouble>     (i)) internal_render(element.getVector().template getAbs<NamedDouble>     (i), child_element);
        else if(element.getVector().template isType<NamedInt>        (i)) internal_render(element.getVector().template getAbs<NamedInt>        (i), child_element);
        else if(element.getVector().template isType<NamedString>     (i)) internal_render(element.getVector().template getAbs<NamedString>     (i), child_element);
        else if(element.getVector().template isType<RenderingRules>  (i)) internal_render(element.getVector().template getAbs<RenderingRules>  (i), child_element);
        else if(element.getVector().template isType<Section>         (i)) internal_render(element.getVector().template getAbs<Section>         (i), child_element);
        else if(element.getVector().template isType<SpannedCell>     (i)) internal_render(element.getVector().template getAbs<SpannedCell>     (i), child_element);
        else if(element.getVector().template isType<TableColumn>     (i)) internal_render(element.getVector().template getAbs<TableColumn>     (i), child_element);
        else if(element.getVector().template isType<TableRow>        (i)) internal_render(element.getVector().template getAbs<TableRow>        (i), child_element);
        else if(element.getVector().template isType<UnavailableCell> (i)) internal_render(element.getVector().template getAbs<UnavailableCell> (i), child_element);
      }
    }
};

} // End namespace OUTPUT
} // End namespace SAGE


#endif
