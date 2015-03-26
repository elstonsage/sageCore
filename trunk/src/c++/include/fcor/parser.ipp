// ---------------------------------------------------------------------------
// Inline Implementation of FcorParser
// ---------------------------------------------------------------------------

inline string
FcorParser::get_weight(bool intraclass) const
{
  weight_type w = my_analysis_options.class_weight;

  if( w == PAIR_WISE )
    return "Equal Weight to Pairs";

  if( w == UNIFORM ) 
    return "Uniform Weight to Pedigrees";

  if( w == MEAN ) 
    return "Mean Weight";

  return "Quadratic Weight Method";
}

inline const RPED::RefMultiPedigree*
FcorParser::get_multi_pedigree() const
{
  return my_multipedigree;
}

inline const vector<name_index_type>&
FcorParser::get_trait_list() const
{
  return my_trait_list;
}

inline const analysis_option_type&
FcorParser::get_analysis_options() const
{
  return my_analysis_options;
}

inline const FcorParser::name_map_type&
FcorParser::get_name_map() const
{
  return my_relationship_name_map;
}

inline size_t
FcorParser::get_trait_count() const
{
  return my_trait_list.size();
}

inline size_t
FcorParser::get_pedigree_count() const
{
  return my_multipedigree->pedigree_count();
}

