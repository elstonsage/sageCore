#ifndef OUTPUT_TABLE_H
#define OUTPUT_TABLE_H

#include <math.h>
#include <string>
#include <sstream>
#include "numerics/functions.h"
#include "error/internal_error.h"
#include "output/Element.h"
#include "output/Misc.h"
#include "output/RenderingRules.h"
#include "output/Graph.h"

namespace SAGE   {
namespace OUTPUT {

/// \internal
/// \brief A spanned cell (only allowed to be created by the TableRow, via its spanLatestCell() function.
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
class SpannedCell : public Element<SpannedCell> 
{
  friend class TableRow;

private:
  SpannedCell() : Element<SpannedCell>("SpannedCell") {}  
};

/// \brief An unavailable cell
///
/// ALLOWABLE CHILD ELEMENTS: [None]
///
/// If a cell is unavailable for any reason, insert an UnavailableCell instance to so indicate.
/// By default, unavailable cells will indicate to the user in some way that they are
/// in fact unavailable values. If, however, you simply want the cell skipped over, instantiate
/// your UnavailableCell with visible = false. Rendering devices will interpret this to mean
/// that the cell should be left blank, rather than use the column's unavailable code.
///
/// \par XML Attributes
///
/// "VISIBLE" (String - "TRUE" / "FALSE") - See above explanation
///
class UnavailableCell : public Element<UnavailableCell, AttrMgr<HasVisible> > 
{
  typedef Element<UnavailableCell, AttrMgr<HasVisible> > Base;

public:

  /// @name Constructors
  //@{
  
    explicit UnavailableCell(bool visible = true) : Base("UnavailableCell")
    {
      setVisible(visible);
    }
    
    UnavailableCell(const UnavailableCell & other) : Base(other) {}
    
    UnavailableCell & operator= (const UnavailableCell & other)
    {
      Base::operator=(other);
      
      return *this;
    }
    
  //@}
};

/// \brief Describes a single row in a table
///
/// ALLOWABLE CHILD ELEMENTS: [SpannedCell, UnavailableCell, String, Double, Int]
///
/// \par Inserting cells
///
/// Please note that if you insert a Double whose value is not-a-number or not finite,
/// the Double will NOT be inserted; instead, the string "Irrational value" will be inserted.
///
/// \par Cell spanning
///
/// It is possible to span a cell across several fields. To do this, you invoke the
/// spanLatestCell() function immediately after having inserted the cell you want to span.
/// This will "pad out" the row by adding SpannedCell instances to fill out the padded
/// fields. The cell content (in pretty print rendering) will be center-justified in the
/// spanned section.
///
/// For instance, let's say you've got a row that must have 5 cells in it (because the 
/// Table in question has five TableColumn's). But you want the middle three cells
/// merged as a single cell. Then you could enter the following code:
///
/// \code
/// TableRow r;
/// r.insert(Double(1.0);
/// r.insert(Int(5);
/// r.spanLatestCell(3);
/// r.insert(String("Hello!"));
/// my_table.insert(r);
/// \endcode
///
class TableRow : public Element<TableRow, AttrMgr<HasGroupName>, SpannedCell, UnavailableCell, String, Double, Int>
{
  typedef Element<TableRow, AttrMgr<HasGroupName>, SpannedCell, UnavailableCell, String, Double, Int> BaseType;

public:

  /// For spanning a cell.

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    TableRow() : BaseType("TableRow") {}

    ///
    /// Copy constructor.
    TableRow(const TableRow & other) : BaseType(other) {}
    
    ///
    /// Assignment operator.
    TableRow operator= (const TableRow & other) { BaseType::operator=(other);  return *this; }
    
  //@}
  
  /// @name Cell-spanning
  //@{
  
    ///
    /// Returns the spanned width of a given cell.
    ///
    int getSpan(size_t column_idx) const
    {
      int width = 1;

      while(1)
      {
        if(column_idx + width >= getVector().size())
          break;

        if(getVector().isType<SpannedCell> (column_idx + width))
        {
          width++;
        }
        else
        {
          break;
        }
      }
        
      return width;      
    }
  
    ///
    /// Span the most recent added cell to a total column width.
    void spanLatestCell(int total_spanned_width)
    {
      for(int i = 0; i < total_spanned_width - 1; ++i)
        insert(SpannedCell());
    }
    
    typedef Command1<int, TableRow, &TableRow::spanLatestCell> SPAN_LATEST_CELL;
    
  //@}
};

/// \brief Represents the presence of a single column in table
///
/// ALLOWABLE CHILD ELEMENTS: [RenderingRules]
///
/// Please note that the table column exists only to indicate
/// that a column is \b present in a table. The actual data
/// is stored in TableRow's which are added directly to the table.
/// Please see the Table reference for more information.
///
/// \par Justification
///
/// For text cells (String's, that is), you can specify the justification for output-rendering
/// purposes. The default justification is LEFT (see Justification for a complete list).
/// The functions for getting / setting justification are getJustification() and setJustification().
///
/// \par Unavailable code
///
/// For any given column there is an unavailable code (a string sequence indicating the value
/// is not available). This will appear if any cell is of type UnavailableCell, and that cell's
/// "visible" option is true (see UnavailableCell for more information). By default, the
/// unavailable code is "Unavailable".
///
/// \par XML Attributes
///
/// "JUSTIFICATION" - (String, "CENTER" / "LEFT" / "RIGHT") - Justification for output
///
class TableColumn : public Element<TableColumn, AttrMgr<HasTitle, HasJustification, HasUnavailableCode, HasGroupName>, RenderingRules>
{
  typedef Element<TableColumn, AttrMgr<HasTitle, HasJustification, HasUnavailableCode, HasGroupName>, RenderingRules> BaseType;

public:

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    /// \param title The name of the column
    /// \param j The justification for any String cells in the column.
    /// \param unavailable_code The string code used for UnavailableCell (if they're visible)
    ///
    explicit TableColumn(const std::string & title, HasJustification::Justification j = LEFT, const std::string & unavailable_code = "Unavailable") : BaseType("TableColumn") 
    { 
      setTitle           (title); 
      setJustification   (j); 
      setUnavailableCode (unavailable_code);
    }

    ///
    /// Copy constructor.
    TableColumn(const TableColumn & other) : BaseType(other) {}
    
    ///
    /// Assignment operator.
    TableColumn& operator= (const TableColumn & other) { BaseType::operator=(other);  return *this; }
    
  //@}
};

/// \brief Describes a set of data in tabular form
///
/// ALLOWABLE CHILD ELEMENTS: [TableColumn, TableRow, RenderingRules, Graph, NamedString]
///
/// \par Introduction
///
/// A Table can hold one or more TableColumn's, TableRow's, and Graph's. You can add them in any
/// order.
///
/// \par TableColumn
///
/// If you want a column to have a title, you should add a TableColumn
/// instance.
///
/// TableColumn's can be grouped by using the Table's beginColumnGroup()
/// and endColumnGroup() functions. When beginColumnGroup() is invoked,
/// you indicate the name of the group in question. All subsequent inserts of
/// TableColumn's will tag those columns as belonging to the group.
///
/// When endColumnGroup() is invoked, subsequent inserts of TableColumn's will
/// tag those columns as not having a group name.
///
/// Keep in mind that if you add a TableRow with more cells in it than there are
/// columns in the table, the table will have the requisite TableColumn's added to it!
///
/// \par TableRow
/// 
/// When you add a TableRow, if it has fewer cells in it than the number of columns
/// already in the Table, the TableRow will have UnavailableCell's added to it to pad
/// out the remaining portion.
///
/// If the TableRow has more cells than columns, the table will have unnamed columns
/// added to it, and all the other rows will be padded out with UnavailableCell's to
/// account for the change in column count.
///
/// TableRow's can be grouped by using the Table's beginRowGroup()
/// and endRowGroup() functions. When beginRowGroup() is invoked,
/// you indicate the name of the group in question. All subsequent inserts of
/// TableRow's will tag those rows as belonging to the group.
///
/// When endRowGroup() is invoked, subsequent inserts of TableRow's will
/// tag those rows as not having a group name.
///
/// \par Table comments
///
/// You can attach an arbitrary number of comments to a Table by inserting
/// NamedString, NamedInt, and NamedDoble instances. 
/// The renderer will generally place all table comments
/// at the end of the rendered table in the order in which you inserted
/// the elements.
///
/// \par Pretty-print rendering logic
///
/// Overall, the pretty print rendering of the Table tries to guess what you have
/// in mind for your table. If your Table's title is a blank string, for instance,
/// no title will be printed for your table. If all of your columns have blank
/// strings for both column group name and column title, then no column headers
/// will be printed.
///
/// \par Graphics rendering
///
/// You can add as many default graph-rendering options as you want to any Table.
/// Each graph is specified by inserting a Graph instance into your table. Please
/// see the Graph reference for more information.
///
/// \par Runtime output
///
/// If you want your table to produce output during runtime whenever content is
/// added, you can do this with the enableRuntimeOutput() function.
///
/// When you invoke this function, you pass it the output stream to which
/// output should be sent, and an (optional) default field width. With
/// pretty-printing, a table is able to precalculate the necessary field widths. With
/// runtime output on-the-fly, however, the table cannot know such information in
/// advance. Consequently, the runtime output will format all fields of equal width
/// (given by the field width argument to enableRuntimeOutput() ). In addition,
/// any comments added to the table will not show up in runtime output, because
/// the program has no way of knowing when the table has been finished.
///
/// You can also disable runtime output with the disableRuntimeOutput() function.
///
/// \par XML Attributes
///
/// The column groupings are indicated by the attribute "COLUMN_GROUPS" which 
/// is a comma-delimited sequence of strings. For every i'th indexed value in
/// the sequence, the string in question indicates the group to which the i'th 
/// column belongs.
///
/// The row groupings are indicated by the attribute "ROW_GROUPS" which 
/// is a comma-delimited sequence of strings. For every i'th indexed value in
/// the sequence, the string in question indicates the group to which the i'th 
/// row belongs.
///
class Table : public Element<Table, AttrMgr<HasTitle, HasTooltip, HasHelpText>, RenderingRules, TableColumn, TableRow, Graph, NamedString, NamedDouble, NamedInt>
{
  typedef Element<Table, AttrMgr<HasTitle, HasTooltip, HasHelpText>, RenderingRules, TableColumn, TableRow, Graph, NamedString, NamedDouble, NamedInt> BaseType;

public:

  friend class ElementBase<Table>;

  /// @name Constructors
  //@{
  
    ///
    /// Constructor.
    /// \param title The title for the table (default is blank).
    explicit Table(const std::string & title = "");

    ///
    /// Copy constructor.
    Table(const Table & other);
    
    ///
    /// Assignment operator.
    Table & operator=(const Table & other);

  //@}
  
  /// @name Column grouping
  //@{
  
    ///
    /// Indicates that all TableColumn's subsequently inserted are to belong
    /// to a column grouping of the given name.
    void beginColumnGroup(const std::string & group_name) { my_current_column_group = group_name; }
    typedef Command1<const std::string &, Table, &Table::beginColumnGroup> BEGIN_COLUMN_GROUP;
    
    ///
    /// Indicates that all TableColumn's subsequently inserted are NOT to belong
    /// to any column grouping.
    void endColumnGroup() { my_current_column_group = ""; }
    typedef Command0<Table, &Table::endColumnGroup> END_COLUMN_GROUP;
  
    ///
    /// Inserts a TableRow composed entirely of UnavailableCell instances (whose 'visible' is
    /// set to \c false).
    void insertBlankRow();
    typedef Command0<Table, &Table::insertBlankRow> INSERT_BLANK_ROW;

    ///
    /// This adds a single-line text message as the next table row. It is equivalent to:
    /// \code
    /// beginRowGroup(message);
    /// insertBlankRow();
    /// endRowGroup();
    /// \endcode
    void insertRowMessage(const std::string & message);
    typedef Command1<const std::string &, Table, &Table::insertRowMessage> INSERT_ROW_MESSAGE;

  //@}

  /// @name Row grouping
  //@{
  
    ///
    /// Indicates that all TableRow's subsequently inserted are to belong
    /// to a column grouping of the given name.
    void beginRowGroup(const std::string & group_name) { my_current_row_group = group_name; }
    typedef Command1<const std::string &, Table, &Table::beginRowGroup> BEGIN_ROW_GROUP;
    
    ///
    /// Indicates that all TableRow's subsequently inserted are NOT to belong
    /// to any column grouping.
    void endRowGroup() { my_current_row_group = ""; }
    typedef Command0<Table, &Table::endRowGroup> END_ROW_GROUP;
  
  //@}

  /// @name Runtime output
  //@{
  
    ///
    /// Tells the table to print out content information to the given std::ostream
    /// as data is added to it.
    ///
    /// \param o The ostream to which data is sent
    /// \param width The width for all fields
    void enableRuntimeOutput(std::ostream * o = &std::cout);
    
    ///
    /// Tells the table to cease printing out content information as it is added.
    void disableRuntimeOutput();
  
  //@}

private:

  ///
  /// Indicates whether or not this Table is supposed to produce runtime output as
  /// it is populated with data.
  bool generateRuntimeOutput() const { return my_ostream != NULL; }

  //=====================================
  // State variables for data population
  //=====================================

  std::string my_current_column_group;
  std::string my_current_row_group;

  //=================================================
  // State variables for output formatting at runtime
  //=================================================

  mutable std::ostream * my_ostream;
  static const size_t s_width = 16;
};

//====================
//  Table validators
//====================

ELEMENT_VALIDATOR(Table, TableColumn)
{
  // Add an entry for the column-group to which this column belongs:

  child.setGroupName(parent.my_current_column_group);
  
  // Add invisible unavailable cells to rows that need them:

  int column_cnt = parent.getVector().count<TableColumn> () + 1;
  
  for(AnyVectorItrs::Iterator<TableRow> row = parent.getVector().begin<TableRow> (); row != parent.getVector().end<TableRow> (); ++row)
  {
    int cells_needed = column_cnt - row->getVector().size();

    for(int i = 0; i < cells_needed; ++i)
      row->insert(UnavailableCell(false));
  }

  // Return success:
  
  return true;
}

ELEMENT_VALIDATOR(Table, Graph)
{
  // Verify the X-Column is available:

  if(!parent.hasElement<TableColumn>(child.getXColumn()))
    return false;
  
  // Verify that every Y-Column is available:
  
  for(Graph::ColumnList::const_iterator i = child.beginYColumn(); i != child.endYColumn(); ++i)
    if(!parent.hasElement<TableColumn>(*i))
      return false;

  return true;
}

ELEMENT_VALIDATOR(Table, TableRow)
{
  // Add an entry for the row-group to which this row belongs:
  
  child.setGroupName(parent.my_current_row_group);

  int cell_cnt   = child.getVector().size(),
      column_cnt = parent.getVector().count<TableColumn>();
      
  if(cell_cnt < column_cnt) // If there aren't enough cells in the row, add invisble UnavailableCell's:
  {
    for(int i = 0; i < column_cnt - cell_cnt; ++i)
      child.insert(UnavailableCell(false));
  }
  else if(cell_cnt > column_cnt) // Add additional columns and cells to all rows:
  {
    for(int i = 0; i < cell_cnt - column_cnt; ++i)
      parent.insert(TableColumn(""));
  }
  
  // Check for runtime output:
  
  if(parent.my_ostream)
  {
    size_t row_cnt = parent.getVector().count<TableRow>() + 1;

    if(row_cnt == 1) // Print out column headers:
    {
      for(AnyVectorItrs::ConstIterator<TableColumn> col_itr  = parent.getVector().begin <TableColumn> ();
                                                    col_itr != parent.getVector().end   <TableColumn> (); ++col_itr)
      {
        if(col_itr->getTitle().length() <= Table::s_width)
        {
          *parent.my_ostream << std::setw(Table::s_width + 1) << std::left << col_itr->getTitle();
        }
        else
        {
          *parent.my_ostream << std::setw(Table::s_width + 1) << std::left << col_itr->getTitle().substr(0, Table::s_width);
        }
      }
      
      *parent.my_ostream << std::endl << std::endl;
    }
    
    for(size_t x = 0; x < child.getVector().size(); ++x)
    {
      std::ostringstream s;
      
      if(child.getVector().isType<String>(x)) 
      {
        s << child.getVector().getAbs<String>(x).toVal();
      }
      else if(child.getVector().isType<Double>(x)) 
      {
        double d = child.getVector().getAbs<Double>(x).toVal();
        
        if(SAGE::isnan(d))
        {
          s << "Non-number";
        }
        else if(finite(d))
        {
          double             d2 = d == 0 ? 0 : d;
          std::ostringstream temp;
          
          temp << d2;
          
          s << (d2 >= 0 ? " " : "") 
            << std::setprecision(temp.str().find("e") == std::string::npos ? 6 : 2) 
            << (temp.str().find("e") == std::string::npos ? std::fixed : std::scientific) 
            << d2;

          // Strip the scientific exponent of any leading zeroes (if necessary):
          s.str(UTIL::StringUtils::stripLeadingExpZero(s.str()));
        }
        else
        {
          s << (d > 0 ? "" : "-") << "Infinity";
        }
      }
      else if(child.getVector().isType<Int>(x)) 
      {
        int d = child.getVector().getAbs<Int>(x).toVal();
        
        if(SAGE::isnan(d))
          s << "Non-number";
        else if(finite(d))
          s << (d >= 0 ? " " : "") << d;
        else
          s << (d > 0 ? "" : "-") << "Infinity";
      }
      
      std::ostringstream s2;
      
      s2 << std::setw(Table::s_width) << std::left << (s.str().length() <= Table::s_width ? s.str() : s.str().substr(0, Table::s_width));
      
      *parent.my_ostream << s2.str() << " ";
    }
    
    *parent.my_ostream << std::endl;

  } // End runtime-output portion

  return true;
}


} // End namespace OUTPUT
} // End namespace SAGE

#endif
