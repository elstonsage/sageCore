#include "output/ViewPlaintext.h"
#include "output/ViewXML.h"
#include "output/ViewPrettyPrint.h"

namespace SAGE   {
namespace OUTPUT {

//===================================
//
//  Section
//
//===================================
std::string 
ViewPrettyPrint::internal_render(const Section & section, const UTIL::OutlineCntr & _cntr)
{    
  // Section title / header:    

  std::ostringstream s;

  if(section.getTitle() != "")
  {
    std::ostringstream hbar;
    std::string        complete_title = section.getTitle();

    hbar << std::setw(10 + complete_title.length()) << std::setfill('=') << "" << std::setfill(' ')  << std::endl;

    s << hbar.str() << std::endl << std::setw(5) << "" << complete_title << std::endl << std::endl << hbar.str() << std::endl;
  }

  UTIL::OutlineCntr cntr(_cntr);
  cntr.increaseDepth();
      
  // Children: Section, Table, String, NamedDouble, NamedInt, NamedString, List
  for(size_t i = 0; i < section.getVector().size(); ++i)
  {
         if(section.getVector().isType<Table>           (i)) s << internal_render(section.getVector().getAbs<Table>       (i), cntr++);    
    else if(section.getVector().isType<String>          (i)) s << internal_render(section.getVector().getAbs<String>      (i), cntr++);    
    else if(section.getVector().isType<List>            (i)) s << internal_render(section.getVector().getAbs<List>        (i), cntr++);    
    else if(section.getVector().isType<NamedDouble>     (i)) s << internal_render(section.getVector().getAbs<NamedDouble> (i), cntr++);    
    else if(section.getVector().isType<NamedInt>        (i)) s << internal_render(section.getVector().getAbs<NamedInt>    (i), cntr++);    
    else if(section.getVector().isType<NamedString>     (i)) s << internal_render(section.getVector().getAbs<NamedString> (i), cntr++);    
    else if(section.getVector().isType<Section>         (i)) s << internal_render(section.getVector().getAbs<Section>     (i), cntr++);    
  }    

  return s.str();    
}     

//=========================================
//
//  List
//
//=========================================
std::string 
ViewPrettyPrint::internal_render(const List & list_, const UTIL::OutlineCntr & _cntr, bool is_top_level_list, HasBulletType::BulletType e)
{
  // Grab and process the outline counter text:
  std::ostringstream        s;
  UTIL::OutlineCntr         cntr;
  std::ostringstream        cntr_txt;
  HasBulletType::BulletType t(list_.getBulletType() == HasBulletType::UNSPECIFIED ? e : list_.getBulletType());
  
  if(!is_top_level_list)
  {
    cntr = _cntr;
  }
    
  if(cntr.getDepth() > 0)
  {
    if(t == HasBulletType::BULLETED)
    {
      cntr_txt << std::setw(cntr.toString(1).length() - 1) << "" << "* ";
    }
    else // t == NUMBERED
    {
      cntr_txt << cntr.toString(1) << " ";
    }
  }

  // If the depth is greater than 0, or if the depth is 0 and the title if present, deliver the text:
  if((cntr.getDepth() > 0) || (cntr.getDepth() == 0 && list_.getTitle() != ""))
  {
    s << cntr_txt.str() 

      << UTIL::StringUtils::indentString(
           UTIL::StringUtils::lineWrapString(
             list_.getTitle(), 
             screen_width - cntr_txt.str().length()), 
           cntr_txt.str().length())
         .substr(cntr_txt.str().length())

      << std::endl 
      << std::endl;
  }
    
  // Increase the depth:
  cntr.increaseDepth();

  // Process children:
  for(size_t i = 0; i < list_.getVector().size(); ++i)
  {
    if(list_.getVector().isType<List>(i))
    {
      s << internal_render(list_.getVector().getAbs<List>(i), cntr++, false, t);
    }
    else // Is String
    {
      s << UTIL::StringUtils::indentString(
           UTIL::StringUtils::lineWrapString(
           list_.getVector().getAbs<String>(i).toVal(), screen_width - cntr_txt.str().length()), cntr_txt.str().length()) 
        << std::endl
        << std::endl;
    }
  }
    
  // Return string:
  return s.str();
}

//===========================================
//
//  Table
//
//===========================================
std::string 
ViewPrettyPrint::internal_render(const Table & table, const UTIL::OutlineCntr & _cntr)
{
  return TableRenderer(table).render(_cntr);
}

//======================================
//
//  TableRenderer::Static variables
//
//======================================
std::string ViewPrettyPrint::TableRenderer::delimiter            = "  ";
std::string ViewPrettyPrint::TableRenderer::horizontal_delimiter = "--";
char        ViewPrettyPrint::TableRenderer::horizontal_char      = '-';

//==========================================
// TableRenderer::CONSTRUCTORS
//==========================================
ViewPrettyPrint::TableRenderer::TableRenderer(const Table         & table) : my_table(table)          { }
ViewPrettyPrint::TableRenderer::TableRenderer(const TableRenderer & other) : my_table(other.my_table) { }

//===========================================
//
// TableRenderer::render(...)
//
//===========================================
std::string
ViewPrettyPrint::TableRenderer::render(const UTIL::OutlineCntr & _cntr) const
{
  // Reset state variables:
  my_left_lengths   . clear();
  my_rendered_cells . clear();
  
  // Calculate formatted field widths and produce output:
  calculateLeftLengths   ();
  calculateRenderedCells ();
  calculateMaxWidths     ();

  std::string s = generateHeader() + generateColumns() + generateRows() + generateFooter();
  
  // Return string:

  return s;
}

//====================================
//
//  calculateLeftLengths()
//
//====================================
void
ViewPrettyPrint::TableRenderer::calculateLeftLengths() const
{
  // Figure out the greatest absolute values for columns with Double's:

  const RenderingRules & table_rules = my_table.getVector().hasOne<RenderingRules>() ? my_table.getVector().getOnly<RenderingRules>() : default_rules;

  my_left_lengths.resize(my_table.getVector().count<TableColumn> (), 0);

  for(size_t column_idx = 0; column_idx < my_left_lengths.size(); ++column_idx)
  {
    const TableColumn    & column       = my_table.getVector().getRel<TableColumn>(column_idx);
    const RenderingRules & column_rules = column.getVector().hasOne<RenderingRules>() ? column.getVector().getOnly<RenderingRules>() : table_rules;

    for(AnyVectorItrs::ConstIterator<TableRow> row = my_table.getVector().begin<TableRow> (); row != my_table.getVector().end<TableRow> (); ++row)
    {
      if(!row->getVector().isType<Double> (column_idx))
        continue;

      if(column_rules.getNumberFormat() == RenderingRules::SCIENTIFIC)
      {
        my_left_lengths[column_idx] = 1; 
        break;
      }

      if(row->getVector().getAbs<Double> (column_idx).getVType() != HasValue<double>::RATIONAL)
        continue;

      size_t l = getLeftLength(row->getVector().getAbs<Double> (column_idx).toVal());

      if(l > my_left_lengths[column_idx])
        my_left_lengths[column_idx] = l;
    }
  }
}

//=======================================
//
//  calculateRenderedCells()
//
//=======================================
void
ViewPrettyPrint::TableRenderer::calculateRenderedCells() const
{
  int column_cnt = my_table.getVector().count<TableColumn>();

  // ColumnIdx x RowIdx:
  my_rendered_cells.resize(column_cnt, std::vector<std::string>(my_table.getVector().count<TableRow>(), ""));
  
  // Fetch table rules:
  const RenderingRules & table_rules = my_table.getVector().hasOne<RenderingRules>() ? my_table.getVector().getOnly<RenderingRules>() : default_rules;

  // Iterate across columns:
  size_t column_idx = 0;

  for(AnyVectorItrs::ConstIterator<TableColumn> column_itr  = my_table.getVector().begin <TableColumn>(); 
                                                column_itr != my_table.getVector().end   <TableColumn>(); ++column_itr, ++column_idx)
  {
    const RenderingRules & column_rules = column_itr->getVector().hasOne<RenderingRules>() ? column_itr->getVector().getOnly<RenderingRules>() : table_rules;
      
    // Iterate across rows:
    size_t row_idx = 0;
    
    for(AnyVectorItrs::ConstIterator<TableRow> row_itr  = my_table.getVector().begin <TableRow>(); 
                                               row_itr != my_table.getVector().end   <TableRow>(); ++row_itr, ++row_idx)
    {
      my_rendered_cells[column_idx][row_idx] = renderCell(*row_itr, column_idx, column_rules);
    }
  }
}

//======================================
//
//  renderCell()
//
//======================================
std::string 
ViewPrettyPrint::TableRenderer::renderCell(const TableRow & row, size_t column_idx, const RenderingRules & rules) const
{
  std::ostringstream s;

  if(row.getVector().isType<String>(column_idx)) 
  {
    s << row.getVector().getAbs<String>(column_idx);
  }
  else if(row.getVector().isType<Double>(column_idx)) 
  {
    const Double & cell = row.getVector().getAbs<Double>(column_idx);
    
    if(cell.getVType() != HasValue<double>::RATIONAL)
    {
      s << rules.render(cell);
    }
    else
    {
      if(!rules.exceedsThreshold(cell.getValue())) // Pad out the left margin
      {
        s << std::setw(my_left_lengths[column_idx] - getLeftLength(cell.getValue())) << std::setfill(' ') << "";
      }

      s << rules.render(cell);
    }
  }
  else if(row.getVector().isType<Int> (column_idx)) 
  {
    s << rules.render(row.getVector().getAbs<Int>(column_idx));
  }
  else if(row.getVector().isType<SpannedCell>(column_idx)) 
  {
    // 
  }
  else if(row.getVector().isType<UnavailableCell>(column_idx)) 
  {
    s << (row.getVector().getAbs<UnavailableCell>(column_idx).getVisible() ? my_table.getVector().getRel<TableColumn>(column_idx).getUnavailableCode() : "");
  }
  else 
  {
    // 
  }

  return s.str();
}

//======================================
//
//  calculateMaxWidths(...)
//
//====================================== 
void 
ViewPrettyPrint::TableRenderer::calculateMaxWidths() const
{
  // Now move on to the columns:

  my_widths.clear();
  my_widths.resize(my_table.getVector().count<TableColumn> ());

  for(int column_idx = my_widths.size() - 1; column_idx >= 0; --column_idx)
  {
    // Find out the widths of the *succeeding* columns (for doing cell spans):

    size_t column_cnt = my_table.getVector().count<TableColumn>();

    // Fetch the column & rendering rules:

    const TableColumn & column = my_table.getVector().getRel<TableColumn>(column_idx);

    // Figure out about the column name and column group name:

    std::string column_title = column.getTitle(),
                group_name   = column.getGroupName();
    size_t      max_width    = column_title.length();

    // If there is a group name and it is longer than the column title:

    if((group_name != "") && (group_name.length() > column_title.length()))
    {
      bool is_first_column        = column_idx == 0,
           is_first_of_this_group = is_first_column || (!is_first_column && (my_table.getVector().getRel<TableColumn>(column_idx - 1).getGroupName() != group_name));

      // If this is the first column of the column group:

      if(is_first_of_this_group)
      {
        // Figure out how wide the entire set of columns in this group is:

        size_t remaining_width = 0;

        for(size_t i = column_idx + 1; i < column_cnt; ++i)
          if(my_table.getVector().getRel<TableColumn>(i).getGroupName() == group_name)
            remaining_width += delimiter.length() + my_widths[i];

        size_t available_width = column_title.length() + remaining_width;

        if(group_name.length() > available_width)
          max_width += group_name.length() - remaining_width;
      }
    }

    // Figure out about the column's data:

    for(size_t row_idx = 0; row_idx < my_table.getVector().count<TableRow>(); ++row_idx)
    {
      const TableRow & row = my_table.getVector().getRel<TableRow>(row_idx);
      
      size_t cell_width = my_rendered_cells[column_idx][row_idx].length(),
             span       = row.getSpan(column_idx);
           
      if(span > 1)
      {
        size_t available_space = max_width;
      
        for(size_t i = 1; i < span; ++i)
          available_space += delimiter.length() + my_widths[column_idx + i];
        
        if(cell_width > available_space)
          max_width += cell_width - available_space;
      }
      else if(cell_width > max_width)
      {
        max_width = cell_width;
      }
    }

    my_widths[column_idx] = max_width;

  } // End of column-loop
}

//================================
//
//  generateHeader()
//
//================================
std::string
ViewPrettyPrint::TableRenderer::generateHeader() const
{
  std::ostringstream s;
  
  if(my_table.getTitle() != "")
  {
    s << std::setw(my_table.getTitle().length() + 10) << std::setfill('=') << "" << std::endl
      << "     "
      << my_table.getTitle() << std::endl
      << std::setw(my_table.getTitle().length() + 10) << std::setfill('=') << "" << std::endl
      << std::setfill(' ')
      << generateBlankLine();
  }
  
  return s.str();
}

//========================================
//
//  generateBlankLine(...)
//
//========================================
std::string
ViewPrettyPrint::TableRenderer::generateBlankLine() const
{
  return "\n";
}

//=======================================
//
//  generateHBar(...)
//
//=======================================
std::string
ViewPrettyPrint::TableRenderer::generateHBar() const
{
  std::ostringstream horizontal_line_sstream;

  for(size_t i = 0; i < my_table.getVector().count<TableColumn>(); ++i)
  {
    horizontal_line_sstream << (!i ? horizontal_delimiter : "")
                            << std::setw(my_widths[i])
                            << std::setfill(horizontal_char) 
                            << "" 
                            << horizontal_delimiter;
  }

  horizontal_line_sstream << std::endl;
 
  return horizontal_line_sstream.str();
}

//========================================
//
//  generateColumns(...)
//
//========================================
std::string
ViewPrettyPrint::TableRenderer::generateColumns() const
{
  std::ostringstream s;
  
  // Column titles / group names:

  size_t column_cnt             = my_table.getVector().count<TableColumn> ();
  bool   has_column_group_names = false,
         has_column_titles      = false;

  for(AnyVectorItrs::ConstIterator<TableColumn> column_itr = my_table.getVector().begin<TableColumn>(); column_itr != my_table.getVector().end<TableColumn>(); ++column_itr)
  {
    if(column_itr->getTitle() != "")
      has_column_titles = true;

    if(column_itr->getGroupName() != "")
      has_column_group_names = true;
  }

  // Anything to report?

  if(!has_column_titles && !has_column_group_names)
    return s.str();
    
  s << generateHBar();

  if(has_column_group_names)
  {
    s << delimiter;

    for(size_t column_idx = 0; column_idx < column_cnt; ++column_idx)
    {
      std::string group_name = my_table.getVector().getRel<TableColumn>(column_idx).getGroupName();
    
      if(group_name == "") // No group name!
      {
        s << std::setw(my_widths[column_idx]) << "" << delimiter;
      }
      else // Group name!
      {
        bool is_first_column   = column_idx == 0,
             is_first_of_group = is_first_column || (!is_first_column && (group_name != my_table.getVector().getRel<TableColumn>(column_idx - 1).getGroupName()));

        if(is_first_of_group)
        {
          size_t field_width = my_widths[column_idx];
      
          for(size_t j = column_idx + 1; j < my_table.getVector().count<TableColumn> (); ++j)
          {
            if(group_name == my_table.getVector().getRel<TableColumn>(j).getGroupName())
            {
              field_width += delimiter.length() + my_widths[j];
            }
          }
          
          size_t pad_width   = field_width - group_name.length(),
                 left_width  = pad_width / 2,
                 right_width = left_width + (pad_width % 2 ? 1 : 0);
        
          s << std::setw(left_width) << "" << group_name << std::setw(right_width) << "" << delimiter;
        }
      }
    }

    s << std::endl;  

  } // End if-has-column-group-names

  if(has_column_titles)
  {
    s << delimiter;

    for(size_t i = 0; i < column_cnt; ++i)
    {
      const std::string & title = my_table.getVector().getRel<TableColumn>(i).getTitle();

      s << std::left << std::setw(my_widths[i]) << title;
      
      s << delimiter;
    }

    s << std::endl;

  } // End if-has-column-titles

  s << generateHBar();

  // Return string:
  return s.str();
}

//======================================
//
//  generateRow(...)
//
//======================================
std::string
ViewPrettyPrint::TableRenderer::generateRow(size_t row_idx, const TableRow & row, const RenderingRules & table_rules) const
{
  // Set up ostringstream:
  
  std::ostringstream s;
  
  // Fetch row:

  // Check row group name:

  bool row_group_is_different = true;
  
  if(row_idx > 0)
  {
    row_group_is_different = row.getGroupName() != my_table.getVector().getRel<TableRow>(row_idx - 1).getGroupName();
  }
  
  if(row_group_is_different)
  {
    if(row_idx)
      s << generateBlankLine();

    if(row.getGroupName() != "")
      s << delimiter << row.getGroupName() << std::endl;
  }

  // Loop through cells and display them:

  s << delimiter;

  int span = 0;

  for(size_t column_idx = 0; column_idx < row.getVector().size(); column_idx += span)
  {
    // Fetch column specific rendering rules:

    const TableColumn    & column       = my_table.getVector().getRel<TableColumn>(column_idx);
    const RenderingRules & column_rules = column.getVector().hasOne<RenderingRules>() ? column.getVector().getOnly<RenderingRules>() : table_rules;

    // Process cell span:

    span = row.getSpan(column_idx);

    int width = -delimiter.length();

    for(size_t i = column_idx; i < column_idx + span; ++i)
      width += my_widths[i] + delimiter.length();

    // Fetch cell text and figure out margin / justification, and deliver it:
    bool        is_spanned                     = span > 1,
                is_double_cell                 = row.getVector().isType<Double> (column_idx),
                double_past_threshold          = is_double_cell && column_rules.exceedsThreshold(row.getVector().getAbs<Double> (column_idx).toVal()),
                treat_as_string_cell           = row.getVector().isType<String> (column_idx) || double_past_threshold,
                unspanned_centered_string_cell = !is_spanned && treat_as_string_cell && (column.getJustification() == TableColumn::CENTER),
                unspanned_left_string_cell     = !is_spanned && treat_as_string_cell && (column.getJustification() == TableColumn::LEFT),
                center_justify                 = is_spanned || unspanned_centered_string_cell,
                left_justify                   = (!is_spanned && is_double_cell) || unspanned_left_string_cell;
    std::string cell_text                      = my_rendered_cells[column_idx][row_idx];

    int remaining_width = width - cell_text.length(),
        left_width      = center_justify ? remaining_width / 2                : (left_justify ? 0 : remaining_width),
        right_width     = center_justify ? left_width + (remaining_width % 2) : (left_justify ? remaining_width : 0);

    s << std::setw(left_width) << std::setfill(' ') << "" << cell_text << std::setw(right_width) << std::setfill(' ') << "" << delimiter;
  }

  s << std::endl;

  return s.str();
}

//=======================================
//
//  generateFooter()
//
//=======================================
std::string
ViewPrettyPrint::TableRenderer::generateFooter() const
{
  std::ostringstream s;

  bool comments_present = false;

  for(size_t i = 0; i < my_table.getVector().size(); ++i)
  {
    // Create a value which is a mini type id of the correct type.  It will
    // be >0 if it's a Named type, 0 otherwise
    size_t type_val = 0;
  
         if(my_table.getVector().isType<NamedDouble> (i)) type_val = 1;
    else if(my_table.getVector().isType<NamedInt>    (i)) type_val = 2;
    else if(my_table.getVector().isType<NamedString> (i)) type_val = 3;

    // If type_val is non-zero, we want to print our comments, if any
    
    if(type_val)
    {
      if(!comments_present)
        s << generateHBar() << generateBlankLine();

      comments_present = true;

      switch(type_val)
      {
        case 1 : /* NamedDouble */ s << my_table.getVector().getAbs<NamedDouble>(i); break;
        case 2 : /* NamedInt    */ s << my_table.getVector().getAbs<NamedInt>   (i); break;
        case 3 : /* NamedString */ s << my_table.getVector().getAbs<NamedString>(i); break;
      }
    }
  }
  
  s << std::endl;

  return s.str();
}

//==================================
//
//  generateRows(...)
//
//==================================
std::string
ViewPrettyPrint::TableRenderer::generateRows() const
{
  std::ostringstream s;
  
  const RenderingRules & table_rules = my_table.getVector().hasOne<RenderingRules>() ? my_table.getVector().getOnly<RenderingRules>() : default_rules;

  for(size_t row_idx = 0; row_idx < my_table.getVector().count<TableRow>(); ++row_idx)
    s << generateRow(row_idx, my_table.getVector().getRel<TableRow>(row_idx), table_rules);
    
  return s.str();
}



} // End namespace OUTPUT
} // End namespace SAGE
