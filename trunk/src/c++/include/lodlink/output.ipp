//============================================================================
// File:      output.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/27/3 - created.                                djb
//                                                                          
// Notes:     inline implementations of classes use by output routines.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  col_hdr
//============================================================================
//
inline
col_hdr::labels::labels(const string& one, const string& two, const string& three)
{
  data.push_back(&one);
  data.push_back(&two);
  data.push_back(&three);
}

inline
col_hdr::col_hdr(size_t label_width, size_t spacer_width, const labels& ls)
      : my_label_width(label_width), my_spacer_width(spacer_width), my_labels(ls)
{
  for(int i = 0; i < 3; ++i)
  {
    assert((my_labels.data[i])->size() <= my_label_width);
  }
}      

inline size_t
col_hdr::lw() const
{
  return my_label_width;
}

inline size_t
col_hdr::sw() const
{
  return my_spacer_width;
}

inline size_t
col_hdr::tw() const
{
  return my_label_width + my_spacer_width;
}

inline const string&
col_hdr::operator[](size_t i) const
{
  assert(i < my_labels.data.size());
  return  *(my_labels.data[i]);
}


//============================================================================
// IMPLEMENTATION:  header
//============================================================================
//
inline
header::header(size_t offset, size_t underline)
      : my_offset(offset), my_underline(underline)
{}

inline size_t
header::offset() const
{
  return  my_offset;
}

// - 'net' width of ith column.
//
inline size_t
header::col_w(size_t i) const
{
  assert(i < my_cols.size());
  return  (my_cols[i])->lw();
}

// - spacer width of ith column.
//
inline size_t
header::spc_w(size_t i) const
{
  assert(i < my_cols.size());
  return  (my_cols[i])->sw();
}

inline bool
header::underline() const
{
  return  my_underline != '0';
}

inline const col_hdr&
header::operator[](size_t index) const
{
  assert(index < my_cols.size());
  return *(my_cols[index]);
}

inline void
header::set_offset(size_t offset)
{
  my_offset = offset;
}

inline void
header::set_underline(char underline)
{
  my_underline = underline;
}

inline void
header::add_col(const col_hdr& col)
{
  my_cols.push_back(&col);
}



