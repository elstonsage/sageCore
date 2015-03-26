//============================================================================
// File:      parser.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/5/4 created        -djb
//                                                                          
// Notes:     Inline implementation of class, parser.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
inline const instructions&
parser::user_instructions() const
{
  return my_instructions;
}

inline size_t
parser::analysis_id() const
{
  return my_analysis_id;
}

inline string
parser::trait_usage_2_string(parser::trait_usage usage)
{
  string  result = "";

  switch(usage)
  {
    case FOR_FAMILY_REP:
      result = "family_rep";
      break;
      
    case FOR_PARTITIONING:
      result = "partition";
      break;
      
    case FOR_POOL_SIZE:
      result = "pool_size_trait";
      break;
      
    case FOR_ALLELE_PCT:
      result = "allele";
      break;
      
    default:
      assert(false);
  }
  
  return  result;
}

inline bool  
parser::trait_registered(size_t t) const
{
  map<size_t, trait_info>::const_iterator  search_result = my_trait_registry.find(t);
  map<size_t, trait_info>::const_iterator  end_iter      = my_trait_registry.end();
  
  return  search_result != end_iter;
}
  
inline bool  
parser::string_registered(size_t s) const
{
  map<size_t, trait_info>::const_iterator  search_result = my_string_registry.find(s);
  map<size_t, trait_info>::const_iterator  end_iter      = my_string_registry.end();
  
  return  search_result != end_iter;
}

inline void  
parser::write_trait_registry(ostream& out) const
{
  map<size_t, trait_info>::const_iterator  iter     = my_trait_registry.begin();
  map<size_t, trait_info>::const_iterator  end_iter = my_trait_registry.end();
  for(; iter != end_iter; ++iter)
  {
    out << "Index " << iter->first << "  ";
    write_trait_info(out, iter->second);
  }
}


//============================================================================
// IMPLEMENTATION:  parser::duplication
//============================================================================
//
inline
parser::duplication::duplication(const pair<string, value>& sub_pop)
      : my_sub_pop(sub_pop)
{}

inline bool
parser::duplication::operator ()(const pair<string, value>& elem)
{
  return  elem.first  == my_sub_pop.first  ||
          elem.second == my_sub_pop.second   ;
}

inline void  
parser::write_trait_info(ostream& out, const parser::trait_info& ti)
{
  out << "trait name " << ti.name << "  trait usage " << trait_usage_2_string(ti.usage) << endl;
}

//============================================================================
// IMPLEMENTATION:  parser::trait_info
//============================================================================
//
inline
parser::trait_info::trait_info(const string& n, parser::trait_usage u)
    : name(n), usage(u)
{}

