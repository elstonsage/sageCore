// ---------------------------------------------------------------------------
// Inline Implementation of SubAVCalT
// ---------------------------------------------------------------------------

inline
size_t
SubAVCalT::pedigree_count() const
{
  return my_pairsetdata->multi_pedigree()->pedigree_count();
}

inline
size_t
SubAVCalT::trait_count() const
{
  return my_pairsetdata->fcor_parser()->trait_list().size();
}

inline
size_t
SubAVCalT::xy_pair_count(size_t ped) const
{
  return (*my_pairset_xy)[ped].size();
}

inline
void
SubAVCalT::pedigree_weight(TriangleMatrix<internal_real_type>& result, size_t p) const
{
  for( size_t i = 0; i < xy_pair_count(p); ++i )
    for( size_t m = 0; m < result.linear_size(); ++m )
      result[m] += my_weight_w[p][i][m];
}

inline
double
SubAVCalT::get_corinfo(size_t i, size_t j, size_t p, size_t t1, size_t t2) const
{
  size_t p1, p2, mem1, mem2;

  p1   = (*my_pairset_xy)[p][i].member_pair.first ->pedigree()->index();
  p2   = (*my_pairset_xy)[p][j].member_pair.second->pedigree()->index();
  mem1 = (*my_pairset_xy)[p][i].member_pair.first ->index();
  mem2 = (*my_pairset_xy)[p][j].member_pair.second->index();

  if( p1 != p2 )
    return std::numeric_limits<double>::quiet_NaN();

  return my_correlation->get_cor_map_cache()[p1](mem1, mem2)->second.correlation(t1, t2);
}

// end of SubAVCalT Implementation





