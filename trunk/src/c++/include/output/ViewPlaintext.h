#ifndef OUTPUT_VIEW_PLAINTEXT_H
#define OUTPUT_VIEW_PLAINTEXT_H

#include "output/Graph.h"
#include "output/Misc.h"
#include "output/RenderingRules.h"
#include "output/Table.h"
#include "util/OutlineCntr.h"

namespace SAGE   {
namespace OUTPUT {

class ViewPlaintext
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

    ///
    /// Plaintext rendering (mostly for debugging purposes).
    template<typename T>
    static std::string render(const T & element, const UTIL::OutlineCntr & _cntr = UTIL::OutlineCntr())
    {
      std::ostringstream s;

      UTIL::OutlineCntr cntr(_cntr);

      // Name:

      s << UTIL::StringUtils::indentString(element.getType(), cntr.getDepth() * 2);

      // Attributes:
 
      typename T::AttributeMap attr_map;

      element.populateAttrMap(attr_map);
      
      if(attr_map.size() > 0)
      {
        s << ", ";

        for(typename T::AttributeMap::const_iterator i = attr_map.begin(); i != attr_map.end(); ++i)
          s << i->first << "=\"" << i->second << "\" ";
      }
      
      s << std::endl;

      // Child elements:

      cntr.increaseDepth();
      
      for(size_t i = 0; i < element.getVector().size(); ++i)
      {
        // Table, Double, Int, String, Graph, List, NamedDouble, NamedInt, NamedString, RenderingRules, Section, SpannedCell, TableColumn, TableRow, UnavailableCell
        
             if(element.getVector().template isType<Table>           (i)) s << ViewPlaintext::render(element.getVector().template getAbs<Table>           (i), cntr++);
        else if(element.getVector().template isType<Double>          (i)) s << ViewPlaintext::render(element.getVector().template getAbs<Double>          (i), cntr++);
        else if(element.getVector().template isType<Int>             (i)) s << ViewPlaintext::render(element.getVector().template getAbs<Int>             (i), cntr++);
        else if(element.getVector().template isType<String>          (i)) s << ViewPlaintext::render(element.getVector().template getAbs<String>          (i), cntr++);
        else if(element.getVector().template isType<Graph>           (i)) s << ViewPlaintext::render(element.getVector().template getAbs<Graph>           (i), cntr++);
        else if(element.getVector().template isType<List>            (i)) s << ViewPlaintext::render(element.getVector().template getAbs<List>            (i), cntr++);
        else if(element.getVector().template isType<NamedDouble>     (i)) s << ViewPlaintext::render(element.getVector().template getAbs<NamedDouble>     (i), cntr++);
        else if(element.getVector().template isType<NamedInt>        (i)) s << ViewPlaintext::render(element.getVector().template getAbs<NamedInt>        (i), cntr++);
        else if(element.getVector().template isType<NamedString>     (i)) s << ViewPlaintext::render(element.getVector().template getAbs<NamedString>     (i), cntr++);
        else if(element.getVector().template isType<RenderingRules>  (i)) s << ViewPlaintext::render(element.getVector().template getAbs<RenderingRules>  (i), cntr++);
        else if(element.getVector().template isType<Section>         (i)) s << ViewPlaintext::render(element.getVector().template getAbs<Section>         (i), cntr++);
        else if(element.getVector().template isType<SpannedCell>     (i)) s << ViewPlaintext::render(element.getVector().template getAbs<SpannedCell>     (i), cntr++);
        else if(element.getVector().template isType<TableColumn>     (i)) s << ViewPlaintext::render(element.getVector().template getAbs<TableColumn>     (i), cntr++);
        else if(element.getVector().template isType<TableRow>        (i)) s << ViewPlaintext::render(element.getVector().template getAbs<TableRow>        (i), cntr++);
        else if(element.getVector().template isType<UnavailableCell> (i)) s << ViewPlaintext::render(element.getVector().template getAbs<UnavailableCell> (i), cntr++);
      }
      
      // Return it all:

      return s.str();
    }
};

} // End namespace OUTPUT
} // End namespace SAGE


#endif
