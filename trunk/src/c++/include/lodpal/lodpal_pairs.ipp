////////////////////////////////////////////////////////////////////////////
//             Implementation of lodpal_pairs (Inline)                    //
////////////////////////////////////////////////////////////////////////////

//
//--------------------------------------------------------------------------
//

inline
pair_filter&
lodpal_pairs::filter()
{
  return my_filter;
}

inline
double lodpal_pairs::member_trait(const MPED::member_base* m, size_t t) const
{
  return pairs.trait(*m, t);
}

inline
bool
lodpal_pairs::built_pairs_info() const
{
  return my_built_pairs_info;
}

inline
bool
lodpal_pairs::built_pairs_info_x() const
{
  return my_built_pairs_info_x;
}

inline
bool
lodpal_pairs::re_built_pairs_info() const
{
  return my_re_built_pairs_info;
}

inline
bool
lodpal_pairs::removed_biggest_pair() const
{
  if(    my_removed_fsib_pairs.size()
      || my_removed_hsib_pairs.size()
      || my_removed_other_pairs.size() )
    return true;
 
  return false;
}

inline
bool
lodpal_pairs::is_mm_pair(size_t prior_index) const
{
  if( !built_pairs_info_x() )
    return false;

  if(    prior_index == 0
      || prior_index == 3
      || prior_index == 6
      || prior_index == 10
      || prior_index == 17 )
    return true;

  return false;
}

inline
bool
lodpal_pairs::is_mf_pair(size_t prior_index) const
{
  if( !built_pairs_info_x() )
    return false;

  if( is_mm_pair(prior_index) || is_ff_pair(prior_index) )
    return false;

  return true;
}

inline
bool
lodpal_pairs::is_ff_pair(size_t prior_index) const
{
  if( !built_pairs_info_x() )
    return false;

  if(    prior_index == 2
      || prior_index == 5
      || prior_index == 9
      || prior_index == 13
      || prior_index == 15
      || prior_index == 16
      || prior_index == 19
      || prior_index == 20 )
    return true;

  return false;
}

inline
bool
lodpal_pairs::is_ff_sib_pair(size_t prior_index) const
{
  if( prior_index == 2 )
    return true;

  return false;
}

inline
bool
lodpal_pairs::parent_of_origin_allowed() const
{
  return my_parent_of_origin_allowed;
}

inline
size_t
lodpal_pairs::re_built_covariate() const
{
  return my_re_built_covariate;
}

inline
size_t
lodpal_pairs::pair_count() const
{
  size_t removed_pair =   my_removed_other_pairs.size()
                        + my_removed_fsib_pairs.size()
                        + my_removed_hsib_pairs.size();

  return my_pair_count - removed_pair;
}

inline
size_t
lodpal_pairs::fsib_pair_count() const
{
  return my_fsib_pair_count - my_removed_fsib_pairs.size();
}

inline
size_t
lodpal_pairs::hsib_pair_count() const
{
  return my_hsib_pair_count - my_removed_hsib_pairs.size();
}

inline
size_t
lodpal_pairs::other_pair_count() const
{
  return pair_count() - fsib_pair_count() - hsib_pair_count();
}

inline
size_t
lodpal_pairs::mm_pair_count() const
{
  size_t removed_mm_pair = 0;

  for( size_t i = 0; i < my_removed_other_pairs.size(); ++i )
    if( is_mm_pair(pairs_info()[my_removed_other_pairs[i]].prior_x_ibd_index) )
      removed_mm_pair += 1;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_mm_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_mm_pair += 1;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_mm_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_mm_pair += 1;

  return my_mm_pair_count - removed_mm_pair;
}

inline
size_t
lodpal_pairs::mm_fsib_pair_count() const
{
  size_t removed_mm_fsib_pair = 0;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_mm_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_mm_fsib_pair += 1;

  return my_mm_fsib_pair_count - removed_mm_fsib_pair;
}

inline
size_t
lodpal_pairs::mfm_hsib_pair_count() const
{
  size_t removed_mfm_hsib_pair = 0;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_mm_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_mfm_hsib_pair += 1;

  return my_mfm_hsib_pair_count - removed_mfm_hsib_pair;
}

inline
size_t
lodpal_pairs::mm_other_pair_count() const
{
  return mm_pair_count() - mm_fsib_pair_count() - mfm_hsib_pair_count();
}

inline
size_t
lodpal_pairs::mm_invalid_pair_count() const
{
  return my_mm_invalid_pair_count;
}

inline
size_t
lodpal_pairs::mf_pair_count() const
{
  size_t removed_mf_pair = 0;

  for( size_t i = 0; i < my_removed_other_pairs.size(); ++i )
    if( is_mf_pair(pairs_info()[my_removed_other_pairs[i]].prior_x_ibd_index) )
      removed_mf_pair += 1;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_mf_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_mf_pair += 1;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_mf_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_mf_pair += 1;

  return my_mf_pair_count - removed_mf_pair;
}

inline
size_t
lodpal_pairs::mf_fsib_pair_count() const
{
  size_t removed_mf_fsib_pair = 0;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_mf_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_mf_fsib_pair += 1;

  return my_mf_fsib_pair_count - removed_mf_fsib_pair;
}

inline
size_t
lodpal_pairs::mff_hsib_pair_count() const
{
  size_t removed_mff_hsib_pair = 0;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_mf_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_mff_hsib_pair += 1;

  return my_mff_hsib_pair_count - removed_mff_hsib_pair;
}

inline
size_t
lodpal_pairs::mf_other_pair_count() const
{
  return mf_pair_count() - mf_fsib_pair_count() - mff_hsib_pair_count();
}

inline
size_t
lodpal_pairs::mf_invalid_pair_count() const
{
  return my_mf_invalid_pair_count;
}

inline
size_t
lodpal_pairs::ff_pair_count() const
{
  size_t removed_ff_pair = 0;

  for( size_t i = 0; i < my_removed_other_pairs.size(); ++i )
    if( is_ff_pair(pairs_info()[my_removed_other_pairs[i]].prior_x_ibd_index) )
      removed_ff_pair += 1;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_ff_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_ff_pair += 1;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_ff_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_ff_pair += 1;

  return my_ff_pair_count - removed_ff_pair;
}

inline
size_t
lodpal_pairs::ff_fsib_pair_count() const
{
  size_t removed_ff_fsib_pair = 0;

  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    if( is_ff_pair(pairs_info()[my_removed_fsib_pairs[i]].prior_x_ibd_index) )
      removed_ff_fsib_pair += 1;

  return my_ff_fsib_pair_count - removed_ff_fsib_pair;
}

inline
size_t
lodpal_pairs::fff_hsib_pair_count() const
{
  size_t removed_fff_hsib_pair = 0;

  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    if( is_ff_pair(pairs_info()[my_removed_hsib_pairs[i]].prior_x_ibd_index) )
      removed_fff_hsib_pair += 1;

  return my_fff_hsib_pair_count - removed_fff_hsib_pair;
}

inline
size_t
lodpal_pairs::ff_other_pair_count() const
{
  return ff_pair_count() - ff_fsib_pair_count() - fff_hsib_pair_count();
}

inline
size_t
lodpal_pairs::ff_invalid_pair_count() const
{
  return my_ff_invalid_pair_count;
}

inline
size_t
lodpal_pairs::max_ped_name_size() const
{
  return my_max_ped_name;
}

inline
size_t
lodpal_pairs::max_ind_name_size() const
{
  return my_max_ind_name;
}

inline
double
lodpal_pairs::get_drp_prob_share(size_t pt, size_t f) const
{
  if( f == 0 )
    return my_drp_ibd_info[pt].first;

  return my_drp_ibd_info[pt].second;
}

inline
double
lodpal_pairs::get_drp_x_prob_share(size_t pt, size_t f) const
{
  if( f == 0 )
    return my_drp_x_ibd_info[pt].first;

  return my_drp_x_ibd_info[pt].second;
}

inline
lodpal_pairs::pair_info_type&
lodpal_pairs::pairs_info()
{
  return my_pairs_info;
}

inline
const lodpal_pairs::pair_info_type&
lodpal_pairs::pairs_info() const
{
  return my_pairs_info;
}

inline
lodpal_pairs::pair_map_type&
lodpal_pairs::pairs_map()
{
  return my_pairs_map;
}

inline
const lodpal_pairs::pair_map_type&
lodpal_pairs::pairs_map() const
{
  return my_pairs_map;
}

inline
vector<size_t>&
lodpal_pairs::removed_fsib_pairs()
{
  return my_removed_fsib_pairs;
}

inline
const vector<size_t>&
lodpal_pairs::removed_fsib_pairs() const
{
  return my_removed_fsib_pairs;
}

inline
vector<size_t>&
lodpal_pairs::removed_hsib_pairs()
{
  return my_removed_hsib_pairs;
}

inline
const vector<size_t>&
lodpal_pairs::removed_hsib_pairs() const
{
  return my_removed_hsib_pairs;
}

inline
vector<size_t>&
lodpal_pairs::removed_other_pairs()
{
  return my_removed_other_pairs;
}

inline
const vector<size_t>&
lodpal_pairs::removed_other_pairs() const
{
  return my_removed_other_pairs;
}

inline
RelativePairs&
lodpal_pairs::relative_pairs()
{
  return pairs;
}

inline
const RelativePairs&
lodpal_pairs::relative_pairs() const
{
  return pairs;
}

inline
const lodpal_parameters&
lodpal_pairs::parameters() const
{
  return params;
}

inline
lodpal_parameters&
lodpal_pairs::parameters()
{
  return params;
}

inline
lodpal_pairs::prior_x_ibd_vector&
lodpal_pairs::prior_x_ibds()
{
  return my_prior_x_ibd_vector;
}

inline
const lodpal_pairs::prior_x_ibd_vector&
lodpal_pairs::prior_x_ibds() const
{
  return my_prior_x_ibd_vector;
}
