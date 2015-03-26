#include "output/Output.h"
#include "tdtex/TransmissionTables.h"

namespace SAGE  {
namespace TDTEX {

//================================================
//
// CONSTRUCTOR
//
//================================================
TransmissionTable::TransmissionTable(size_t n)
{
  resize(n);
}

//================================================
//
// COPY CONSTRUCTOR   
//
//================================================
TransmissionTable::TransmissionTable(const TransmissionTable & other) :
  my_counts(other.my_counts)
{ }
        
//================================================
//
// TransmissionTable operator=
//
//================================================
TransmissionTable& 
TransmissionTable::operator=(const TransmissionTable & other)
{
  if(this != &other)
  {
    my_counts = other.my_counts;
  }
  
  return *this;
}
            

//================================================
//
// resize(...)
//
//================================================
void
TransmissionTable::resize(size_t n)
{
  my_counts.resize_fill(n, n, 0);
}

//================================================
//
// max_width(...)
//
//================================================
size_t
TransmissionTable::max_width(size_t min_width) const
{
  size_t w = min_width;

  for(size_t i = 0; i < get_counts().rows(); ++i)
  {
    w = max(w, heading(i).size());

    for(size_t j = 0; j < get_counts().rows(); ++j)
    {
      double n = get_counts()(i, j);

      if(n)
        w = max<size_t>(w, (size_t)ceil(log10(n)));
    }
  }

  return w;
}

struct ReorderComparitor
{
  ReorderComparitor() { }
  
  explicit ReorderComparitor(const TransmissionTable * p) : my_p(p) { }
  
  ReorderComparitor(const ReorderComparitor & other) : my_p(other.my_p) { }
  
  ReorderComparitor& operator= (const ReorderComparitor & other)
  {
    if(this != &other)
    {
      my_p = other.my_p;
    }
    
    return *this;
  }

  inline bool operator() (size_t i, size_t j) const
  {
    return my_p->heading(i) < my_p->heading(j);
  }

  const TransmissionTable * my_p;
};

//==================================================
//
//  populate_reordering_vector(...)
//
//==================================================
void 
TransmissionTable::populate_reordering_vector(ReorderingVector & v) const
{
  // Set up local variables:
  typedef std::set<size_t, ReorderComparitor> AlphabeticallySortedSet;

  ReorderComparitor       comp                      (this);
  AlphabeticallySortedSet alphabetically_sorted_ids (comp);

  // Sort the id's alphabetically:
  for(size_t i = 0; i < get_counts().rows(); ++i)
    alphabetically_sorted_ids.insert(i);  

  // Resize the reordering vector:
  v.resize(get_counts().rows());
  
  // For each element in the reordering vector, locate it in the alphabetically sorted set:
  for(size_t i = 0; i < v.size(); ++i)
  {
    size_t found_index = 0;

    for(AlphabeticallySortedSet::const_iterator itr = alphabetically_sorted_ids.begin(); itr != alphabetically_sorted_ids.end(); ++itr, ++found_index)
    {
      if(i == *itr)
      {
        v[i] = found_index;
        continue;
      }
    }
  }
}

//==================================================
//
//  generate_reordered_matrix()
//
//==================================================
TransmissionTable::ReorderedMatrix 
TransmissionTable::generate_reordered_matrix() const
{
  // Set up local variables:
  ReorderedMatrix  m;
  ReorderingVector v;

  // Copy over old (unordered matrix) for the moment:
  m.matrix = get_counts();
  
  // Figure out the reordering:
  populate_reordering_vector(v);

  // Populate the headings vector:
  m.headings.resize(m.matrix.rows());
  
  for(size_t orig_idx = 0; orig_idx < v.size(); ++orig_idx)
  {
    size_t new_idx = v[orig_idx];
    
    m.headings[new_idx] = heading(orig_idx);
  }

  // Transfer the reordered values into the correct locations:
  for(size_t orig_row_idx = 0; orig_row_idx < m.matrix.rows(); ++orig_row_idx)
  {
    size_t new_row_idx = v[orig_row_idx];
    
    for(size_t orig_col_idx = 0; orig_col_idx < m.matrix.cols(); ++orig_col_idx)
    {
      size_t new_col_idx = v[orig_col_idx];
      
      m.matrix(new_row_idx, new_col_idx) = get_counts()(orig_row_idx, orig_col_idx);
    }
  }

  // Return the ReorderedMatrix:
  return m;
}
    

//==================================================
//
//  dump()
//
//==================================================
void
TransmissionTable::dump() const
{
  OUTPUT::Table t("TransmissionTable: " + units() + " dump");
  
  t << OUTPUT::TableColumn("Not-transmitted")
    << OUTPUT::Table::BEGIN_COLUMN_GROUP("Transmitted");

  for(size_t i = 0; i < get_counts().cols(); ++i)
    t << OUTPUT::TableColumn(heading(i));

  t << OUTPUT::Table::END_ROW_GROUP();
    
  for(size_t i = 0; i < get_counts().rows(); ++i)
  {
    OUTPUT::TableRow r;
    
    r << heading(i);
    
    for(size_t j = 0; j < get_counts().cols(); ++j)
    {
      r << get_counts()(i, j);
    }
    
    t << r;
  }  
  
  std::cout << t << std::flush;
}


} // End namespace TDTEX
} // End namespace SAGE
