#ifndef OUTPUT_GRAPH_H
#define OUTPUT_GRAPH_H

#include <math.h>
#include "error/internal_error.h"
#include "output/Element.h"
#include "output/Misc.h"
#include "output/RenderingRules.h"
#include <string>
#include <sstream>

namespace SAGE   {
namespace OUTPUT {

/// \brief Describes a default graph to generate based on a Table's contents.
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
/// \par Introduction
///
/// A graph describes default graphical rendering options to be associated with a Table.
///
/// When you construct your Graph instance, you pass the constructor the following information:
///
/// - The name of the graph
///
/// - The type to generate (line graph, scatter plot, or bar chart)
///
/// - The name of the column that contains the x-values
///
/// - The name of a column that contains y-values (you can subsequently add additional y-columns)
///
/// \par Attributes summary
///
/// "Title" (String) - Name of graph to generate
///
/// "Type" (String - "LINE_GRAPH", "SCATTER_PLOT", "BAR_CHART") - Type of graph to generate
///
/// "XColumn" (String) - Name of column that has x-coordinate values
///
/// "XColumnType" (String - "NUMERICAL", "CATEGORICAL") - Type of values in the x-column
///
/// "YColumns" (comma-delimited strings) - Names of columns that have y-coordinate values
///
class Graph : public Element<Graph, AttrMgr<HasTitle, HasGraphType, HasXColumn, HasXColumnType, HasYColumns> >
{
public:

  typedef Element<Graph, AttrMgr<HasTitle, HasGraphType, HasXColumn, HasXColumnType, HasYColumns> > Base;

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    /// \param title The title of the graph
    /// \param t The type of graph to draw (see Type)
    /// \param x_column The name of the x-column
    /// \param x_type The type of values in the x-column (see XColumnType)
    /// \param y_column The name of the y-column
    Graph(const std::string                 & title, 
                HasGraphType::GraphType       t, 
          const std::string                 & x_column, 
                HasXColumnType::XColumnType   x_type, 
          const std::string                 & y_column) : Base("Graph")
    { 
      setTitle       (title);
      setGraphType   (t);
      setXColumnType (x_type);
      setXColumn     (x_column);
      addYColumn     (y_column);
    }
    
    ///
    /// Copy constructor.
    Graph(const Graph & other) : Base(other) { }
    
    ///
    /// Assignment operator.
    Graph & operator= (const Graph & other) { Base::operator=(other); return *this; }
  
  //@}
};


} // End namespace OUTPUT
} // End namespace SAGE

#endif
