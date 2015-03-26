#include "output/Table.h"

namespace SAGE {
namespace OUTPUT {

//======================================
//
//  Static variables
//
//======================================

const int screen_width = 72;

//======================================
//
//  insertBlankRow()
//
//======================================
void 
Table::insertBlankRow()
{
  TableRow r;
  for(size_t i = 0; i < getVector().count<TableColumn> (); ++i)
    r << UnavailableCell(false);
  insert(r);
}

//======================================
//
//  enableRuntimeOutput(...)
//
//======================================
void 
Table::enableRuntimeOutput(std::ostream * o)
{
  my_ostream = o;
}

//======================================
//
//  disableRuntimeOutput(...)
//
//======================================
void 
Table::disableRuntimeOutput()
{
  my_ostream = NULL;
}

//======================================
//
//  insertRowMessage()
//
//======================================
void 
Table::insertRowMessage(const std::string & message)
{
  beginRowGroup  (message);
  insertBlankRow ();
  endRowGroup    ();
}

//======================================
//
//  CONSTRUCTOR
//
//====================================== 
Table::Table(const std::string & title) : BaseType("Table")
{
  my_current_column_group = "";
  my_current_row_group    = "";
  my_ostream              = NULL;

  setTitle(title);
}

//======================================
//
//  COPY CONSTRUCTOR
//
//====================================== 
Table::Table(const Table & other) : BaseType(other) 
{
  my_current_column_group = other.my_current_column_group;
  my_current_row_group    = other.my_current_row_group;
  my_ostream              = other.my_ostream;
}

//======================================
//
//  ASSIGNMENT OPERATOR
//
//====================================== 
Table& Table::operator=(const Table & other)
{
  BaseType::operator=(other);

  my_current_column_group = other.my_current_column_group;
  my_current_row_group    = other.my_current_row_group;
  my_ostream              = other.my_ostream;

  return *this;
}



} // End namespace OUTPUT
} // End namespace SAGE
